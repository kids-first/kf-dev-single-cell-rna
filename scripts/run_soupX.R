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

# ElbowPlot can provide an estimation of dimensions to use
print(ElbowPlot(seurat_obj))


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

# To visualize the distribution of expression (relative to what would be expected were each cell pure background) across all cells in the data set. 
# When no geneset is provided, the function will try and guess which genes might be useful.
print(plotMarkerDistribution(sc))

# Correcting expression profile
# This will contain a corrected matrix to be used in place of the original table of counts in downstream analyses.
# use roundToInt option to make sure we output integer matrix.
out <- adjustCounts(sc, roundToInt = TRUE) 

# Optional but we can keep original matrix
# Attention as if we do, then for the next step we need to ensure we are using the correct assay, i.e., `RNA_SoupX`
seurat_obj[["RNA_SoupX"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts)
seurat_obj@assays$RNA@counts <- out

# Setup assay for seurat object
DefaultAssay(seurat_obj) <- "RNA_SoupX"
seurat_obj


#########################################################################################
# Investigating changes in expression
# Before proceeding let's have a look at what this has done. 
# We can get a sense for what has been the most strongly decreased 
# by looking at the fraction of cells that were non-zero now set to zero after correction.
# Rename to toc or it gives NULL
toc = filtered_matrix
toc = sc$toc

cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
print(mostZeroed)


# List of gene markers showing on this list are highly specific markers of one cell type or group of cells.
# This is important to note as it may lead to erroneous inferences of potential cell specific genes.
# For example, presence of mitochondrial genes MT-ND4, MT-ND4L or immune cells.

# If on the other hand we focus on genes for which there is a quantitative difference,
print(tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20))
# Then we might notice different gene markers associated with other pathways and cell types.


# Visualizing expression distribution
# Way back at the start, we did a quick visualization to look at how the ratio of IGKC expression to pure soup was distributed. 
# Now that we've corrected our data, we can see how that compares to our corrected data. 
# The function plotChangeMap can help us with this. 
# By default it plots the fraction of expression in each cell that has been deemed to be soup and removed.
print(plotChangeMap(sc, out, "MT-ND4L")) #mtDNA
print(plotChangeMap(sc, out, "CD300E")) #monocytes
print(plotChangeMap(sc, out, "MT-ND4L")) #monocytes
#########################################################################################

# Integrating with downstream tools
# Of course, the next thing you'll want to do is to load this corrected expression matrix into some downstream analysis tool and further analyse the data.

# The corrected matrix can then be used for any downstream analysis in place of the uncorrected raw matrix. 
# If you are using 10X data and would like to save these final counts out in the same format, 
# you can use the DropletUtils write10xCounts function like this

colnames(out) = gsub('-1$', '', colnames(out))
colnames(out) = paste(sample_name, colnames(out), sep=":")

saveRDS(out, paste0(sample_name, ".rds"))
# saveRDS(out, file.path(results_dir, "seurat_object_RNA_SoupX.rds"))

# If we are not using the below output, let's remove it or comment it out
DropletUtils:::write10xCounts(sample_name, out)


