#!/usr/bin/env Rscript
# QC Script that follows https://github.com/hbctraining/scRNA-seq_online/blob/scRNAseq/lessons/04_SC_quality_control.md

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(scales))
suppressMessages(library(cowplot))
suppressMessages(library(RCurl))
suppressMessages(library(DropletUtils))
suppressMessages(library(optparse))

option_list <- list(
  make_option(
    opt_str = "--data_dir",
    type = "character",
    help = "Path to 10X counts directory, like 'outs/raw_feature_bc_matrix'. Only use if not using h5_counts input"
  ),
  make_option(
    opt_str = "--h5_counts",
    type = "character",
    help = "Path to 10X  h5-formatted counts. Only use if not using matrix counts directory data_dir"
  ),
  make_option(
    opt_str = "--sample_id",
    type = "character",
    help = "Sample name"
  ),
  make_option(
    opt_str = "--output_basename",
    type = "character",
    default = "seurat_qc.filtered",
    help = "Name for filtered Seurat RDS output"
  ),
  make_option(
    opt_str = "--min_features",
    type = "numeric",
    default = 100,
    help = "minimum number of genes that need to be detected per cell. this is an absolute floor"
  ),
  make_option(
    opt_str = "--min_umi",
    type = "numeric",
    default = 500,
    help = "minimum number of umi for cell-level filtering"
  ),
  make_option(
    opt_str = "--min_genes",
    type = "numeric",
    default = 250,
    help = "minimum number of genes for cell-level filtering"
  ),
  make_option(
    opt_str = "--min_complexity",
    type = "numeric",
    default = 0.8,
    help = "minimum novelty score (log10GenesPerUMI)"
  ),
  make_option(
    opt_str = "--max_mito_ratio",
    type = "numeric",
    default = 0.2,
    help = "Maximum ratio of mitochondrial genes per cell"
  ),
  make_option(
    opt_str = "--min_gene_prevalence",
    type = "numeric",
    default = 10,
    help = "Minimum number of cells a gene must be expressed in to keep after filtering"
  )
)


calculate_metrics <- function(seurat_obj){
  message("Calculating novelty score")
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

  message("Calculating mito ratio")
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
  return(seurat_obj)
}

rename_cols <- function(metadata){
  message("Rename orig.ident -> seq_source, nCount_RNA -> nUMI, nFeature_RNA -> nGene for ease of use")
  metadata <- metadata %>%
          dplyr::rename(seq_source = orig.ident,
                        nUMI = nCount_RNA,
                        nGene = nFeature_RNA)
  return(metadata)
}

plot_qc <- function(out_filename, m1, m2, min_umi, min_genes, min_complexity, max_mito_ratio){
  pdf(out_filename, onefile = TRUE)
  m1$sample <- paste0(m1$sample, "_pre-filter")
  m2$sample <- paste0(m2$sample, "_post-filter")
  metadata = rbind(m1, m2)
  message("Plot Number of transcripts per cell")
  # need to assign then "print" each fig so that it gets output to same pdf
  a <- metadata %>% 
      ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
      geom_density(alpha = 0.2) + 
      scale_x_log10() + 
      theme_classic() +
      ggtitle("Number of Transcripts per Cell") +
      xlab("Number of UMIs") +
      ylab("Cell density") +
      geom_vline(xintercept = min_umi)
  print(a)
  message("Plot number of genes per cell")
  b <- metadata %>% 
      ggplot(aes(color=sample, x=nGene, fill= sample)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() +
      ggtitle("Number of Genes per Cell") +
      xlab("Number of Genes") +
      scale_x_log10() + 
      geom_vline(xintercept = min_genes)
  print(b)
  message("Plot cell complexity using novelty score")
  c <- metadata %>%
      ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
      ggtitle("Novelty Score") +
      geom_density(alpha = 0.2) +
      theme_classic() +
      geom_vline(xintercept = min_complexity)
  print(c)
  message("Plot mitochondrial counts ratio")
  d <- metadata %>% 
      ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
      ggtitle("Mito Counts Ratio") +
      geom_density(alpha = 0.2) +
      scale_x_log10() +
      theme_classic() +
      geom_vline(xintercept = max_mito_ratio)
  print(d)
  message("Plot joint filtering effects")
  e <- metadata %>% 
      ggplot(aes(x=nUMI, y=nGene, color=mitoRatio, group=interaction(nUMI, nGene))) + 
      geom_point() + 
      scale_colour_gradient(low = "gray90", high = "black") +
      stat_smooth(method=lm) +
      scale_x_log10() +
      ggtitle("Joint Filter Effects") +
      xlab("Number UMIs") +
      scale_y_log10() +
      ylab("Number of Genes") +
      theme_classic() +
      geom_vline(xintercept = min_umi) +
      geom_hline(yintercept = min_genes) +
      facet_wrap(~sample)
  print(e)
  dev.off()
}

boxplot_summary <- function(m1, m2, table_fn){
  summary_box_stats <- boxplot(m1$nUMI, m2$nUMI, m1$nGene, m2$nGene, m1$log10GenesPerUMI, m2$log10GenesPerUMI, m1$mitoRatio, m2$mitoRatio,
      names=c("Num_UMIs_pre-filter", "Num_UMIs_post-filter", "Num_Genes_pre-filter", "Num_Genes_post-filter", "Novelty_Score_pre-filter", "Novelty_Score_post-filter","Mito_Ratio_pre-filter", "Mito_Ratio_post-filter"))
  col_names = c("Lower_whisker", "Lower_hinge", "Median", "Upper_hinge", "Upper_whisker")
  summary_to_print <- t(summary_box_stats$stats)
  colnames(summary_to_print) <- col_names
  rownames(summary_to_print) <- summary_box_stats$names
  # Do some formatting because R hates you
  formatted_summary <- rownames_to_column(round(as.data.frame(summary_to_print),digits=2), var="Metrics")
  write.table(formatted_summary, table_fn, quote=FALSE, row.names=FALSE, sep="\t")
}

opts <- parse_args(OptionParser(option_list = option_list), print_help_and_exit = TRUE)

message("Loading 10X counts matrix and creating Seurat object")
if (!is.null(opts$data_dir)){
  counts <- Read10X(data.dir = opts$data_dir)
  pname = basename(opts$data_dir)
} else {
  counts <- Read10X_h5(opts$h5_counts, use.names = TRUE)
  pname = basename(opts$h5_counts)
}
seurat_counts <- CreateSeuratObject(counts = counts,
                          min.features = opts$min_features, project=pname)

seurat_counts <- calculate_metrics(seurat_counts)

message("Collating metadata")
metadata_prefiltered <- seurat_counts@meta.data
metadata_prefiltered <- rename_cols(metadata_prefiltered)
metadata_prefiltered$sample <- opts$sample_id
message("Printing prefiltered QC Metrics")
write.table(rownames_to_column(metadata_prefiltered, var="cell_ids"), paste0(opts$output_basename, ".barcode_qc.metrics.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
seurat_counts@meta.data <- metadata_prefiltered

message("Applying cell-level filters")
filtered_seurat <- subset(x = seurat_counts, 
                         subset= (nUMI >= opts$min_umi) & 
                           (nGene >= opts$min_genes) & 
                           (log10GenesPerUMI > opts$min_complexity) & 
                           (mitoRatio < opts$max_mito_ratio))
message("Dropping genes with 0 counts from cell-filtered data")
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= opts$min_gene_prevalence
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
# drop extra columns metadata made by CreateSeuratObject
filtered_seurat@meta.data <- subset(filtered_seurat@meta.data, select = -c(`orig.ident`, `nCount_RNA`, `nFeature_RNA`))

message("Output QC filtered count matrix")
filtered_ct_matrix_fname = paste0(opts$output_basename, ".qc_filtered.counts_matrix.h5")
DropletUtils::write10xCounts(path = filtered_ct_matrix_fname, x = filtered_seurat@assays$RNA@data, barcodes=colnames(filtered_seurat@assays$RNA@data), gene.id=rownames(filtered_seurat@assays$RNA@data), gene.type = "Gene Expression", version=3, type="HDF5", chemistry = "Single Cell 3' v3")
message("Repeating metrics and plots for filtered data")
metadata_clean <- filtered_seurat@meta.data
message("Normalizing read counts")
filtered_seurat <- NormalizeData(filtered_seurat)

# Identify the most variable genes
filtered_seurat <- FindVariableFeatures(filtered_seurat, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
filtered_seurat <- ScaleData(filtered_seurat)
# Identify the 20 most highly variable genes
ranked_variable_genes <- VariableFeatures(filtered_seurat)
top_genes <- ranked_variable_genes[1:20]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 20
message("Plotting variable genes")
p <- VariableFeaturePlot(filtered_seurat)
p <- LabelPoints(plot = p, points = top_genes, repel = TRUE)
var_features_fn <- paste0(opts$output_basename, ".variable_features.pdf")
pdf(var_features_fn, onefile = TRUE)
print(p)
dev.off()

message("Creating QC plots")
plots_fname = paste0(opts$output_basename, ".QC_plots.pdf")
plot_qc(plots_fname, metadata_prefiltered, metadata_clean, opts$min_umi, opts$min_genes, opts$min_complexity, opts$max_mito_ratio)

message("Creating boxplot summary table")
table_fn = paste0(opts$output_basename, ".boxplot_summary_metrics.tsv")
boxplot_summary(metadata_prefiltered, metadata_clean, table_fn)

message("Print before and after cell count")
# calculate individual contributions of cell dropping
umi_low = length(which(seurat_counts$nUMI < opts$min_umi))
gene_low = length(which(seurat_counts$nGene < opts$min_genes))
complexity_low = length(which(seurat_counts$log10GenesPerUMI <= opts$min_complexity))
mito_high = length(which(seurat_counts$mitoRatio >= opts$max_mito_ratio))
ct_head = c("Pre-filter cell counts",
  paste0("UMIs below ", opts$min_umi),
  paste0("Genes below ", opts$min_genes),
  paste0("Complexity below ", opts$min_complexity),
  paste0("Mito ratio above ", opts$max_mito_ratio),
  "Post-filter cell counts")

cell_cts = c(length(seurat_counts$sample), umi_low, gene_low, complexity_low, mito_high, length(filtered_seurat$sample))
ct_tbl_fn <- paste0(opts$output_basename, ".cell_counts.tsv")
ct_tbl = as.data.frame(t(cell_cts))
colnames(ct_tbl) <- ct_head
write.table(ct_tbl, ct_tbl_fn, quote=FALSE, row.names=FALSE, sep="\t")
message("QC complete!")