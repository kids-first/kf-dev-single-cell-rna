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

#process inputs
option_list <- list(
  make_option(
    opt_str = "--raw",
    type = "character",
    help = "Path to h5 file containing raw matrix from Cell Ranger count"
  ),
  make_option(
    opt_str = "--fil",
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
  )
)

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
raw <- opts$raw
fil <- opts$fil
clusters_file <- opts$cluster
sample_name <- opts$sample_name

# load 10X data and convert to Soup Channel (object SoupX uses for analysis)
filtered_matrix <- Read10X_h5(fil, use.names = TRUE)
raw_matrix <- Read10X_h5(raw, use.names = TRUE)

############################################################################
# If cluster file available
clusters <- read.csv(clusters_file)
mDat <- data.frame(clusters=clusters$Cluster, row.names=clusters$Barcode)
sc <- SoupChannel(raw_matrix, filtered_matrix, mDat)

############################################################################
# If cluster file NOT available
# Make a Seurat object from the sparce matrix
seurat_obj  <- CreateSeuratObject(counts = filtered_matrix)
seurat_obj 

# Make a “SoupChannel”, the object needed to run SoupX
sc  <- SoupChannel(raw_matrix, filtered_matrix)
print(sc)

# Let's estimate clusters
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() 
#seurat_obj <- SCTransform(seurat_obj, verbose = F)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = TRUE)

# After clustering is obtained, it can be added to the channel using setClusters. 
# setDR is useful for visualizations.
meta <- seurat_obj@meta.data
umap <- seurat_obj@reductions$umap@cell.embeddings
sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
sc  <- setDR(sc, umap)
head(meta)

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




# clean the data
out_matrix <- adjustCounts(sc, roundToInt = T)
colnames(out_matrix) = gsub('-1$', '', colnames(out_matrix))
colnames(out_matrix) = paste(sample_name, colnames(out_matrix), sep=":")

saveRDS(out_matrix, paste0(sample_name, ".rds"))

DropletUtils:::write10xCounts(sample_name, out_matrix)