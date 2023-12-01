# R script to run seurat analysis cluster. Main aim to to feed SoupX clusters
# leverages code from Alex S. seurat_analysis.R script

suppressMessages(library(optparse))
suppressMessages(library('Seurat'))
suppressMessages(library('tools'))


option_list <- list(

  make_option(
    opt_str = "--name",
    default = "test",
    type = "character",
    help = "Project name"
  ),
  make_option(
    opt_str = "--data",
    default = file.path(getwd(), "data.rds"),
    type = "character",
    help = "Input data file"
  ),

  make_option(
    opt_str = "--num_pcs",
    default = 10,
    type = "numeric",
    help = "Minimum number of principal components to retain for clustering"
  ),
  make_option(
    opt_str = "--pc_cut",
    default = 0.05,
    type = "numeric",
    help = "p-value cutoff for determing the number of prinicipal components to retain for clustering"
  ),
  make_option(
    opt_str = "--knn_granularity",
    default = 0.5,
    type = "numeric",
    help = "KNN clustering granularity parameter"
  ),
  make_option(
    opt_str = "--min_cells",
    default = 3,
    type = "numeric",
    help = "Minimum number of cells a feature must be present in to be retained"
  )
)

# get opts
opts <- parse_args(OptionParser(option_list = option_list))
data_file <- file.path(opts$data)
project_name <- opts$name
norm_method <- opts$norm_method
num_pcs <- opts$num_pcs
pc_cut <- opts$pc_cut
knn_granularity <- opts$knn_granularity
min_cells <- opts$min_cells
# flexible file loader
ext <- file_ext(data_file)
if (ext == "gz" | ext == "tgz") {
  print("using tar file")
  untar(tarfile = data_file, exdir = "./input_matrix")
  analysis.data <- Read10X(data.dir = paste("./input_matrix",
    list.files("./input_matrix")[1], sep = "/"))
  analysis <- CreateSeuratObject(counts = analysis.data, project = project_name)
} else if (ext == "h5") {
  print("using h5 file")
  analysis.data <- Read10X_h5(data_file)
  analysis <- CreateSeuratObject(counts = analysis.data, project = project_name)
} else if (ext == "rds" | ext == "RDS") {
  print("using rds file")
  analysis <- readRDS(file = data_file)
} else {
  stop("Unrecognized input file type")
}
# PCA prep
analysis <- NormalizeData(analysis, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() 

# run PCA
print("running PCA")
analysis <- RunPCA(analysis, verbose=FALSE)
print("done with PCA")

## run jackstraw to score pcs
analysis <- JackStraw(analysis, num.replicate = 100)
analysis <- ScoreJackStraw(analysis, dims = 1:20)

## determine number of pcs to use for clustering
scores <- JS(object = analysis[["pca"]], slot = "overall")
auto_pcs <- length(which(scores[, "Score"] <= pc_cut))
if (num_pcs < auto_pcs) {
  num_pcs <- auto_pcs
}
## clustering
analysis <- FindNeighbors(analysis, dims = 1:num_pcs)
analysis <- FindClusters(analysis, resolution = knn_granularity)

file <- 'clusters.csv'
header <- 'Barcodes,Cluster\n'
cat(header, file = file)
write.table(cbind(rownames(analysis@meta.data),analysis@meta.data$seurat_clusters), file, append=TRUE, sep=',', col.names=FALSE, row.names=FALSE, quote = FALSE)