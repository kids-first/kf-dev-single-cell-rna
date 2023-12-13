# R script for mRNA decontamination of 10X single-cell data using SoupX

options(error = traceback)

#install packages if not installed.
list.of.packages <- c("optparse", "codetools", "survival", "plyr", "zoo", "data.table", "htmlwidgets", "lazyeval", "reshape2", "SoupX")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(new.packages)
if(length(new.packages)) suppressMessages(install.packages(new.packages, repos='http://cran.us.r-project.org'))

suppressMessages(library(optparse))
suppressMessages(library(SoupX))
suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

#process inputs
option_list <- list(
  make_option(
    opt_str = "--raw",
    type = "character",
    help = "Path to h5 file containing raw matrix from Cell Ranger count"
  ),
  make_option(
    opt_str = "--filtered",
    type = "character",
    help = "Path to h5 file containing filtered matrix from Cell Ranger count"
  ),
  make_option(
    opt_str = "--cluster",
    type = "character",
    help = "Path to csv file containing barcodes and cluster id from Cell Ranger count"
  ),
  make_option(
    opt_str = "--sample_name",
    type = "character",
    help = "Name associated with this sample"
  ),
  make_option(
    opt_str = "--results_dir",
    type = "character",
    help = "output dir to put results",
    default="./"
  )
)

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
raw <- opts$raw
filtered <- opts$filtered
clusters_file <- opts$cluster
sample_name <- opts$sample_name
results_dir <- opts$results_dir
dir.create(results_dir)

# load 10X data and convert to Soup Channel (object SoupX uses for analysis)
filtered_matrix <- Read10X_h5(filtered, use.names = TRUE)
raw_matrix <- Read10X_h5(raw, use.names = TRUE)

############################################################################
# SoupX will use default pdf settings, set output file name to be more useful
pdf(file = file.path(results_dir, paste0(sample_name, ".soupx.plots.pdf")))

# Make a Seurat object from the sparce matrix
seurat_obj  <- CreateSeuratObject(counts = filtered_matrix)

# If cluster file available
if (!is.null(clusters_file)) {
  print("Clusters file exists!")
  clusters <- read.csv(clusters_file)
  mDat <- data.frame(clusters=clusters$Cluster, row.names=clusters$Barcode)
  sc <- SoupChannel(raw_matrix, filtered_matrix, mDat)
} else {
  print("Clusters file does not exist. Let's calculate clusters by using Seurat.")  
  # Make a “SoupChannel”, the object needed to run SoupX
  sc  <- SoupChannel(raw_matrix, filtered_matrix)
  print(sc)
  
  # Normalize, find variable features and scale the data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() 
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    
  ## Run jackstraw to score pcs
  seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
  seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
  
  # Define paarmeters for clustering
  pc_cut = 0.05
  num_pcs = 10
  res_value = 0.5
  
  # Determine number of pcs to use for clustering
  scores <- JS(object = seurat_obj[["pca"]], slot = "overall")
  auto_pcs <- length(which(scores[, "Score"] <= pc_cut))
  if (num_pcs < auto_pcs) {
    num_pcs <- auto_pcs
  }
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:num_pcs, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:num_pcs, verbose = FALSE)
  # Estimate clusters
  seurat_obj <- FindClusters(seurat_obj, resolution = res_value, verbose = TRUE)
  
  # ElbowPlot can provide an estimation of dimensions to use
  # print(ElbowPlot(seurat_obj))
  
  # After clustering is obtained, it can be added to the channel using setClusters. 
  # setDR is useful for visualizations.
  meta <- seurat_obj@meta.data
  umap <- seurat_obj@reductions$umap@cell.embeddings
  sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- setDR(sc, umap)
  head(meta)
}


######################################################################################
######################################################################################
# With defined clusters, run the main SoupX function, calculating ambient RNA profile.
######################################################################################
######################################################################################
# to set the contamination fraction to 20% for all cells.
# we shouldn't remove more than 20% background
# sc = setContaminationFraction(sc, 0.2)

# Estimate contamination fraction rho
# The posterior distribution is calculated using a Poisson likelihood with a gamma distribution prior, parametrised by its mean priorRho and standard deviation priorRhoStdDev. The dotted line in the above plot shows the prior distribution. The default parameters have been calibrated to be fairly non-specific 
# with a slight preference towards values of rho in the 0% to 10% range which is most commonly seen for fresh (i.e. not nuclear) single cell experiments.
# The default values place only a very weak constraint, as can be seen by setting a uniform prior.
sc <- autoEstCont(sc) #priorRhoStdDev = 0.3
print(sc)

# Genes with highest expression in background. These are often enriched for ribosomal proteins.
head(sc$metaData)
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)

# To visualize the distribution of expression (relative to what would be expected were each cell pure background) across all cells in the data set. 
# When no geneset is provided, the function will try and guess which genes might be useful.
print(plotMarkerDistribution(sc))

# Save the plot
ggsave(filename = paste0(sample_name, ".soupx.plotMarkerDistribution.pdf"),
       path = results_dir, 
       width = 6, 
       height = 5, 
       device = "pdf", 
       useDingbats = FALSE)

# Correcting expression profile
# This will contain a corrected matrix to be used in place of the original table of counts in downstream analyses.
# use roundToInt option to make sure we output integer matrix.
out <- adjustCounts(sc, roundToInt = TRUE) 

# Let's save corrected matrix as a separate assay
seurat_obj[["RNA_SoupX"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts)
seurat_obj@assays$RNA@counts <- out

# Setup assay for seurat object for next steps to the one with corrected gene counts
DefaultAssay(seurat_obj) <- "RNA_SoupX"


# Integrating with downstream tools
# Of course, the next thing you'll want to do is to load this corrected expression matrix into some downstream analysis tool and further analyse the data.

# The corrected matrix can then be used for any downstream analysis in place of the uncorrected raw matrix. 
# If you are using 10X data and would like to save these final counts out in the same format, 
# you can use the DropletUtils write10xCounts function like this

colnames(out) = gsub('-1$', '', colnames(out))
colnames(out) = paste(sample_name, colnames(out), sep=":")

saveRDS(out, file.path(results_dir, paste0(sample_name, ".soupx.rds")))
