# R script convert Cell Ranger-style counts dir to h5 file output
# leverages code from Alex S. seurat_analysis.R script

suppressMessages(library(optparse))
suppressMessages(library(DropletUtils))

#process inputs
option_list <- list(
  make_option(
    opt_str = "--counts_dir",
    type = "character",
    help = "Path to dir containing raw matrix from Cell Ranger count"
  ),
  make_option(
    opt_str = "--output_basename",
    type = "character",
    help = "Output file basename, like sample.raw"
  ),
  make_option(
    opt_str = "--sample_name",
    type = "character",
    help = "Name associated with this sample"
  )
)

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
counts_dir <- opts$counts_dir
output_basename <- opts$output_basename
sample_name <- opts$sample_name

counts_in = DropletUtils::read10xCounts(counts_dir, sample.names=sample_name, col.names=TRUE)
outfile = paste0(output_basename, ".h5")
DropletUtils::write10xCounts(outfile, counts_in@assays@data@listData$counts, gene.id=counts_in@rowRanges@elementMetadata@listData$ID, gene.symbol=counts_in@rowRanges@elementMetadata@listData$Symbol, type="HDF5")