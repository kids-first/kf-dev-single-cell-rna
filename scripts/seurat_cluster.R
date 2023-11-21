# R script to run seurat analysis cluster. Main aim to to feed SoupX clusters
# leverages code from Alex S. seurat_analysis.R script

suppressMessages(library(optparse))
suppressMessages(library(Seurat))

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
    opt_str = "--out",
    default = file.path(getwd(), "out"),
    type = "character",
    help = "Output directory path"
  ),
  make_option(
    opt_str = "--min_features",
    default = 200,
    type = "numeric",
    help = "Minimum number of genes observed in a cell to retain"
  ),
  make_option(
    opt_str = "--norm_method",
    default = "LogNormalize",
    type = "character",
    help = "Normalization to apply to counts (LogNormalize, CLR, RC)"
  ),
  make_option(
    opt_str = "--retain_features",
    default = 2000,
    type = "numeric",
    help = "Number of most-variable features to initially retain"
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
min_features <- opts$min_features
norm_method <- opts$norm_method
retain_features <- opts$retain_features
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
  analysis <- CreateSeuratObject(counts = analysis.data, project = project_name,
    min.cells = min_cells, min.features = min_features)
} else if (ext == "h5") {
  print("using h5 file")
  analysis.data <- Read10X_h5(data_file)
  analysis <- CreateSeuratObject(counts = analysis.data, project = project_name,
    min.cells = min_cells, min.features = min_features)
} else if (ext == "rds" | ext == "RDS") {
  print("using rds file")
  analysis <- readRDS(file = data_file)
} else {
  stop("Unrecognized input file type")
}
# PCA prep
analysis <- NormalizeData(analysis, normalization.method = norm_method,
scale.factor = 10000)

## identify highly variable genes
analysis <- FindVariableFeatures(analysis, selection.method = "vst",
nfeatures = retain_features)
## scale data
all.genes <- rownames(analysis)
analysis <- ScaleData(analysis, features = all.genes)

# run PCA
print("running PCA")
analysis <- RunPCA(analysis, features = VariableFeatures(object = analysis))
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

header <- 'Barcodes,Cluster'
cat(header, file = file)
write.table(cbind(rownames(analysis@meta.data),analysis@meta.data$seurat_clusters), file, append=TRUE, sep=',', col.names=FALSE, row.names=FALSE, quote = FALSE)