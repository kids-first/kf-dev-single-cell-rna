suppressPackageStartupMessages({
  library(scDblFinder)
  library(Seurat)
  library(scater)
  library(future)
  library(tidyverse)
  library(grid)
  library(optparse)
})

#process inputs
option_list <- list(
  make_option(
    opt_str = "--results_dir",
    type = "character",
    help = "RDS object produced by Seurat QC"
  ),
  make_option(
    opt_str = "--data_path",
    type = "character",
    help = "RDS object produced by SoupX"
  ),
  make_option(
    opt_str = "--sample_name",
    type = "character",
    help = "CSV file with doublet information produced by Scrublet"
  )
)

opts <- parse_args(OptionParser(option_list = option_list), print_help_and_exit = TRUE)

results_dir = opts$results_dir
data_path = opts$data_path
sample_name = opts$sample_name

plan("multisession", workers = 4)
options(future.globals.maxSize = 64000 * 1024^2)

results_dir <-
  file.path(results_dir)
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


# File path to doublets directory
doublets_dir <- file.path(paste0(results_dir,"/scDblFinder-", sample_name, "/"))
if (!dir.exists(doublets_dir)) {
  dir.create(doublets_dir)
}

knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE, 
                      message = FALSE, 
                      fig.path = file.path(results_dir))

seurat_obj_raw <- readRDS(data_path)
DefaultAssay(seurat_obj_raw) <- "RNA"

# Make this reproducible as UMAP algorithm is not deterministic
set.seed(1234)
# Normalize, find variable features and scale data
seurat_obj <- seurat_obj_raw %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() 

# Run dim_reductions
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs=30) %>%
                    RunUMAP(dims = 1:20, n.neighbors=20L, reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)

seurat_obj <- as.SingleCellExperiment(seurat_obj)

# Save sce
saveRDS(seurat_obj, file = paste0(results_dir, "/", "seurat_obj_sce.rds"))

# Run scDblFinder
bp <- BiocParallel::MulticoreParam(4, RNGseed=1234)
BiocParallel::bpstart(bp)
set.seed(1234)
seurat_obj  <- scDblFinder(seurat_obj, clusters = TRUE, BPPARAM = bp)
BiocParallel::bpstop(bp)

# Estimate singlets and doublets in the library
singlets_doublets_number <- print(table(seurat_obj$scDblFinder.class))

# Estimate pct of doublets in the library
doublets_pct_library <- table(seurat_obj$scDblFinder.class)[2]/(table(seurat_obj$scDblFinder.class)[1]+table(seurat_obj$scDblFinder.class)[2])*100

oublets <- seurat_obj@colData[,str_detect(colnames(seurat_obj@colData), "scDblFinder")] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("cell")

# Save table with all predictions
write.csv(doublets, paste0(doublets_dir, "doublets_", sample_name, ".csv"))
    
# Cells that we need to filter from the final object
doublets.to.filter <- subset(doublets, doublets$scDblFinder.class == "doublet")

doublets_number <- print(table(seurat_obj$scDblFinder.class))

# Save table with cells to be filtered
write.csv(doublets.to.filter, paste0(doublets_dir, "doublets_to_filter_", sample_name,".csv"))

# Save seurat_obj with doublets predictions in the metadata
saveRDS(seurat_obj, paste0(doublets_dir, "sce_obj_", sample_name,".rds") )

print("Doublets predictions")
  
fname <- paste0(doublets_dir, sample_name, "-Doublets_prediction.pdf")
print(fname)
print(gridExtra::grid.arrange(
  plotUMAP(object_list[[i]], colour_by="scDblFinder.cluster", point_size = 0.1) + ggtitle("scDblFinder clusters"),
  plotUMAP(object_list[[i]], colour_by="scDblFinder.score", point_size = 0.1) + ggtitle("Doublets score"),
  plotUMAP(object_list[[i]], colour_by="scDblFinder.class", point_size = 0.1) + ggtitle("Doublets prediction"),
  ncol=2,
  top = paste0(sample_name, "-Doublets prediction")))

pdf(file = fname, width = 12, height = 12)

print(gridExtra::grid.arrange(
  plotUMAP(object_list[[i]], colour_by="scDblFinder.cluster", point_size = 0.1) + ggtitle("scDblFinder clusters"),
  plotUMAP(object_list[[i]], colour_by="scDblFinder.score", point_size = 0.1) + ggtitle("Doublets score"),
  plotUMAP(object_list[[i]], colour_by="scDblFinder.class", point_size = 0.1) + ggtitle("Doublets prediction"),
  ncol=2,
  top = paste0(sample_name, "-Doublets prediction")))

dev.off()