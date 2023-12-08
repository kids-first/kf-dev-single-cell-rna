#!/usr/bin/env Rscript

if(!suppressWarnings(require(optparse))){
    install.packages("optparse", quiet=TRUE)
    library(optparse)
}
library(Seurat)

#process inputs
option_list <- list(
  make_option(
    opt_str = "--seurat_raw_rds",
    type = "character",
    help = "Raw RDS object produced by Seurat QC"
  ),
  make_option(
    opt_str = "--seurat_qc_rds",
    type = "character",
    help = "qc RDS object produced by Seurat QC"
  ),
  make_option(
    opt_str = "--soupx_rds",
    type = "character",
    help = "RDS object produced by SoupX"
  ),
  make_option(
    opt_str = "--scdblfinder_tsv",
    type = "character",
    help = "TSV file with doublet information produced by scDblFinder"
  ),
  make_option(
    opt_str = "--output_basename",
    type = "character",
    default = "seurat_qc.filtered",
    help = "Name for filtered Seurat RDS output"
  )
)

opts <- parse_args(OptionParser(option_list = option_list), print_help_and_exit = TRUE)
if (is.null(opts$seurat_qc_rds) | is.null(opts$seurat_qc_rds) | is.null(opts$soupx_rds)) {
    write("\n\nseurat_qc_rds, soupx_rds, and scdblfinder_tsv must be provided! One or more are missing!", stderr())
    quit(status=1)
}
output_basename <- opts$output_basename
seurat_raw_rds <- opts$seurat_raw_rds
seurat_qc_rds <- opts$seurat_qc_rds
soupx_rds <- opts$soupx_rds
scdblfinder_tsv <- opts$scdblfinder_tsv

write("\n\nLoading Files...", stderr())
seurat_raw_obj <- readRDS(seurat_raw_rds)
seurat_qc_obj <- readRDS(seurat_qc_rds)
soupx_obj <- readRDS(soupx_rds)
scDblFinder_df <- read.csv(scdblfinder_tsv, sep = "\t", header = TRUE, row.names=1)
write("Done loading files, getting aligned cell count", stderr())
# Create and update a filter summary file
filter_summary = paste0(output_basename, ".stats.tsv")
write("Filter step\tNumber of to remove\tNumber of cells", file=filter_summary)

num_align <- ncol(seurat_raw_obj)
write(paste("Alignment", "NA", num_align, sep="\t"), file=filter_summary, append = TRUE)

write("Filtering Seurat Cells...", stderr())
write("...Beginning QC Filtration...", stderr())
num_qc_filtered = ncol(seurat_qc_obj)
write(paste("Seurat QC", (num_align - num_qc_filtered), num_qc_filtered, sep="\t"), file=filter_summary, append = TRUE)
# Report mRNA
filter_cells <- colnames(seurat_qc_obj)[which(!(colnames(seurat_qc_obj) %in% colnames(soupx_obj)))] # Get Seurat cells not found in SoupX
if (num_qc_filtered == 0) {
    write("......No Seurat cells identified as mRNA contaminated. No cells removed...", stderr())
} else {
    write("......Cells found in Seurat were identified by SoupX as mRNA contaminated! Removing the following cells from Seurat:", stderr())
    cat(filter_cells, sep = "\n", file = stderr())
}

write("...Beginning SoupX Filtration...", stderr())
# Remove Cells Not Found in Both Seurat and SoupX
common_cells <- Reduce(intersect, list(colnames(seurat_qc_obj),colnames(soupx_obj)))
seurat_qc_obj <- subset(seurat_qc_obj, cells = common_cells)
soupx_obj <- soupx_obj[,common_cells]
# write(paste("SoupX contamination", ncol(seurat_qc_obj), sep="\t"), file=filter_summary, append = TRUE)
write("...SoupX Filtration Complete...", stderr())

# Create new RNA_filter assay
write("...Adding RNA_filter assay to Seurat...", stderr())
seurat_qc_obj[["RNA_filter"]] <- CreateAssayObject(counts = soupx_obj)
DefaultAssay(seurat_qc_obj) <- "RNA_filter"  # make new filter the default for future steps
write("...RNA_filter assay added...", stderr())

# Remove scDblFinder-Identified Doublets from Seurat
write("...Beginning scDblFinder Filtration...", stderr())
seurat_qc_obj <- AddMetaData(seurat_qc_obj, scDblFinder_df)
# Report Doublets
doublet_cells <- colnames(seurat_qc_obj)[which(seurat_qc_obj@meta.data$is_doublet == TRUE)]
dbl_ct = length(doublet_cells)
if (dbl_ct == 0) {
    write("......No Seurat cells identified as doublets. No cells removed...", stderr())
} else {
    write("......Cells found in Seurat were identified by scDblFinder as doublets! Removing the following cells from Seurat:", stderr())
    cat(doublet_cells, sep = "\n", file = stderr())
}

seurat_filtered_obj <- subset(seurat_qc_obj, cells = WhichCells(seurat_qc_obj, expression = is_doublet==FALSE))
write("...scDblFinder Filtration Complete...", stderr())
write("Seurat Filtering Complete!", stderr())

write(paste("scDblFinder doublet", dbl_ct, ncol(seurat_filtered_obj), sep="\t"), file=filter_summary, append = TRUE)

write("Saving Filtered Seurat Object...", stderr())
saveRDS(seurat_filtered_obj, paste0(output_basename, ".rds"))
write("Program Complete!", stderr())
