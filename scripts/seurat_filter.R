#!/usr/bin/env Rscript

suppressMessages(library(optparse))
suppressMessages(library(Seurat))

#process inputs
option_list <- list(
  make_option(
    opt_str = "--seurat_qc_rds",
    type = "character",
    help = "RDS object produced by Seurat QC"
  ),
  make_option(
    opt_str = "--soupx_rds",
    type = "character",
    help = "RDS object produced by SoupX"
  ),
  make_option(
    opt_str = "--scrublet_csv",
    type = "character",
    help = "CSV file with doublet information produced by Scrublet"
  ),
  make_option(
    opt_str = "--output_name",
    type = "character",
    default = "seurat_qc.filtered.rds",
    help = "Name for filtered Seurat RDS output"
  )
)

opts <- parse_args(OptionParser(option_list = option_list), print_help_and_exit = TRUE)
if (is.null(opts$seurat_qc_rds) | is.null(opts$seurat_qc_rds) | is.null(opts$soupx_rds)) {
    write("seurat_qc_rds, soupx_rds, and scrublet_csv must be provided! One or more are missing!", stderr())
    quit(status=1)
}
output_name <- opts$output_name
seurat_qc_rds <- opts$seurat_qc_rds
soupx_rds <- opts$soupx_rds
scrublet_csv <- opts$scrublet_csv

write("Loading Files...", stderr())
seurat_qc_obj <- readRDS(seurat_qc_rds)
soupx_obj <- readRDS(soupx_rds)
scrublet_df <- read.csv(scrublet_csv, sep = ",", header = FALSE, row.names=1)
write("Done loading files!", stderr())

write("Filtering Seurat Cells...", stderr())
write("...Begining SoupX Filtration...", stderr())
# Report mRNA
filter_cells <- colnames(seurat_qc_obj)[which(!(colnames(seurat_qc_obj) %in% colnames(soupx_obj)))] # Get Seurat cells not found in SoupX
if (length(filter_cells) == 0) {
    write("......No Seurat cells identified as mRNA contaminated. No cells removed...", stderr())
} else {
    write("......Cells found in Seurat were identified by SoupX as mRNA contaminated! Removing the following cells from Seurat:", stderr())
    cat(filter_cells, sep = "\n", file = stderr())
}
# Remove Cells Not Found in Both Seurat and SoupX
common_cells <- Reduce(intersect, list(colnames(seurat_qc_obj),colnames(soupx_obj)))
seurat_qc_obj <- subset(seurat_qc_obj, cells = common_cells)
soupx_obj <- soupx_obj[,common_cells]
write("...SoupX Filtration Complete...", stderr())

# Create new RNA_filter assay
write("...Adding RNA_filter assay to Seurat...", stderr())
seurat_qc_obj[["RNA_filter"]] <- CreateAssayObject(counts = soupx_obj)
DefaultAssay(seurat_qc_obj) <- "RNA_filter"  # make new filter the default for future steps
write("...RNA_filter assay added...", stderr())

# Remove Scrublet-Identified Doublets from Seurat
write("...Beginning Scrublet Filtration...", stderr())
colnames(scrublet_df) <- c("doublet_score","is_doublet") # set scrublet colnames
seurat_qc_obj <- AddMetaData(seurat_qc_obj, scrublet_df)

# Report Doublets
doublet_cells <- colnames(seurat_qc_obj)[which(seurat_qc_obj@meta.data$is_doublet == "True")]
if (length(doublet_cells) == 0) {
    write("......No Seurat cells identified as doublets. No cells removed...", stderr())
} else {
    write("......Cells found in Seurat were identified by Scrublet as doublets! Removing the following cells from Seurat:", stderr())
    cat(doublet_cells, sep = "\n", file = stderr())
}

seurat_filtered_obj <- subset(seurat_qc_obj, cells = WhichCells(seurat_qc_obj, expression = is_doublet=="False"))
write("...Scrublet Filtration Complete...", stderr())
write("Seurat Filtering Complete!", stderr())

write("Saving Filtered Seurat Object...", stderr())
saveRDS(seurat_filtered_obj, output_name)
write("Program Complete!", stderr())
