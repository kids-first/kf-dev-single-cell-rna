#!/usr/bin/env Rscript
# Author:  M. Brown
# Date: 2024-06-17
# Description:  Convert h5 files to 10x format for downstream processing with SWANS

suppressMessages(library(optparse))
suppressMessages(library(DropletUtils))
suppressMessages(library(Seurat))

read_h5_output_cr_matrix <- function(h5_file, sample, cr_type) {
  # read in h5 files and write to 10x format
  h5_obj <- Read10X_h5(h5_file, use.names = TRUE)
  cr_matrix_out_path <- file.path(sample, "outs", paste0(cr_type, "_feature_bc_matrix"))
  # Must at least create this path as write10xCounts seems to have a depth limit creating output path
  dir.create(cr_matrix_out_path, recursive = TRUE, showWarnings = FALSE)
  message("Write ", cr_type, " output to ", cr_matrix_out_path)
  write10xCounts(path=cr_matrix_out_path, x=h5_obj, type='auto', version='3', overwrite = TRUE)
}

option_list <- list(
  make_option("--sample",
    type="character",
    help="Sample name as str"
  ),
  make_option("--raw",
    type="character",
    help="raw h5 matrix file from cellranger"
  ),
  make_option("--filtered",
    type="character",
    help="filtered h5 matrix file from cellranger"
  )
)

opt <- parse_args(OptionParser(option_list=option_list))
sample <- if (is.null(opt$sample)) stop("--sample is required. See --help for all opts") else opt$sample
raw_h5 <- if (is.null(opt$raw)) stop("--raw is required. See --help for all opts") else opt$raw
filtered_h5 <- if (is.null(opt$filtered)) stop("--filtered is required. See --help for all opts") else opt$filtered

# read in h5 files and write to 10x format
read_h5_output_cr_matrix(raw_h5, sample, "raw")
read_h5_output_cr_matrix(filtered_h5, sample, "filtered")
