#R script to run seurat analysis

#this should be more manual, but will let us compare
#the expected data to the actual

suppressMessages(library(optparse))
suppressMessages(library("tools"))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))

save_plot <- function(cmd, name) {
  #function to take in a command and a basename and output a plot file
  plot_file <- file.path(out_dir, paste0(name, ".png"))
  png(filename = plot_file, width = plot_size, height = plot_size)
  eval(parse(text = cmd))
}

#process inputs
option_list <- list(
  make_option(
    opt_str = "--size",
    default = 500,
    type = "numeric",
    help = "Size for plots",
  ),
  make_option(
    opt_str = "--name",
    default = "test",
    type = "character",
    help = "Project name"
  ),
  make_option(
    opt_str = "--out_size",
    default = 20,
    type = "numeric",
    help = "Number of genes to include in the cluster output file"
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
    opt_str = "--max_features",
    default = 2500,
    type = "numeric",
    help = "Maximum number of genes observed in a cell to retain"
  ),
  make_option(
    opt_str = "--max_mt",
    default = 5,
    type = "numeric",
    help = "Maximum mitochondrial percentage observed in a cell to retain"
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
    opt_str = "--nheatmap",
    default = 10,
    type = "numeric",
    help = "Number of principal components for which to produce heatmaps"
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
  ),
  make_option(
    opt_str = "--min_samples",
    default = 100,
    type = "numeric",
    help = "Minimum number of samples to continue analysis, used after QC filtering"
  )
)

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
plot_size <- opts$size
clust_size <- opts$out_size
data_file <- file.path(opts$data)
out_dir <- file.path(opts$out)
project_name <- opts$name
min_features <- opts$min_features
max_features <- opts$max_features
max_mt <- opts$max_mt
norm_method <- opts$norm_method
retain_features <- opts$retain_features
nheatmap <- opts$nheatmap
num_pcs <- opts$num_pcs
pc_cut <- opts$pc_cut
knn_granularity <- opts$knn_granularity
min_cells <- opts$min_cells
min_samples <- opts$min_samples

#make output directory
dir.create(out_dir, recursive = "true")

##load data
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

##preprocess / qc
#calculate % MT reads
analysis[["percent.mt"]] <- PercentageFeatureSet(analysis, pattern = "^MT-")

#Visualize QC metrics as a violin plot
name <- "qc_violin"
cmd <- "VlnPlot(analysis, features = c(\"nFeature_RNA\", \"nCount_RNA\",
  \"percent.mt\"), ncol = 3)"
save_plot(cmd, name)

#subset data with desired options
analysis <- subset(analysis, subset = nFeature_RNA > min_features &
  nFeature_RNA < max_features & percent.mt < max_mt)

# output message with qc status if there aren't enough samples
qc_status_file = file.path(out_dir, paste0("qc_status", ".txt"))
if (dim(analysis)[2] < min_samples) {
  # create file with error message
  write("Too few samples remaining after filtering.", file = qc_status_file)
} else {
  write("PASS", file = qc_status_file)

  #normalize data with selected type and scale factor
  analysis <- NormalizeData(analysis, normalization.method = norm_method,
    scale.factor = 10000)

  #identify highly variable genes
  analysis <- FindVariableFeatures(analysis, selection.method = "vst",
    nfeatures = retain_features)

  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(analysis), 10)

  #plot variable features with labels
  plot1 <- VariableFeaturePlot(analysis)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  name <- "features"
  cmd <- "plot2"
  save_plot(cmd, name)

  ##pca

  #scale data
  all.genes <- rownames(analysis)
  analysis <- ScaleData(analysis, features = all.genes)

  #run PCA
  print("running PCA")
  analysis <- RunPCA(analysis, features = VariableFeatures(object = analysis))
  print("done with PCA")

  #create a heat map for the first 10 PCs
  name <- "heat_map"
  cmd <- "DimHeatmap(analysis, dims = 1:nheatmap, cells = 500, balanced = TRUE)"
  save_plot(cmd, name)

  #create an elbow plot for PCA
  name <- "elbow"
  cmd <- "ElbowPlot(analysis)"
  save_plot(cmd, name)

  #run jackstraw to score pcs
  analysis <- JackStraw(analysis, num.replicate = 100)
  analysis <- ScoreJackStraw(analysis, dims = 1:20)

  #determine number of pcs to use for clustering
  scores <- JS(object = analysis[["pca"]], slot = "overall")
  auto_pcs <- length(which(scores[, "Score"] <= pc_cut))
  if (num_pcs < auto_pcs) {
    num_pcs <- auto_pcs
  }

  #output pca summary
  file <- file.path(out_dir, paste0("pca_summary", ".txt"))
  sink(file)
  print(paste("Number of PCs used in clustering:", num_pcs))
  print("PCA Summary")
  print(analysis[["pca"]], dims = 1:7, nfeatures = clust_size)
  sink()

  ##clustering
  analysis <- FindNeighbors(analysis, dims = 1:num_pcs)
  analysis <- FindClusters(analysis, resolution = knn_granularity)

  #run UMAP
  analysis <- RunUMAP(analysis, dims = 1:num_pcs)

  #plot UMAP
  name <- "umap"
  cmd <- "DimPlot(analysis, reduction = \"umap\")"
  save_plot(cmd, name)

  #could also run tSNE at this point

  ##differential expression
  #find markers for every cluster compared to all remaining cells,
  #report only the positive ones
  analysis.markers <- FindAllMarkers(analysis, only.pos = TRUE, min.pct = 0.25,
    logfc.threshold = 0.25)
  file <- file.path(out_dir, paste0("cluster_markers", ".txt"))
  cluster_markers <- analysis.markers %>%
    group_by(cluster) %>%
    top_n(n = clust_size, wt = avg_log2FC)
  write.csv(cluster_markers, file = file)

  #the next few steps are just examples
  #eventually, we'll have to figure out what we would want

  #generate a heat map of the top 10 markers
  top10 <- analysis.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  name <- "cluster_heat"
  cmd <- "DoHeatmap(analysis, features = top10$gene) + NoLegend()"
  save_plot(cmd, name)

  ##save data object
  out_file <- file.path(out_dir, paste0("analysis", ".rds"))
  print("saving final data object")
  saveRDS(analysis, file = out_file)
}
