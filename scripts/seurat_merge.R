# R script for merging 10X single-cell datasets using Seurat

options(error = traceback)

suppressMessages(library(optparse))
suppressMessages(library(Seurat))

#process inputs
option_list <- list(
  make_option(
    opt_str = "--matrix_files",
    type = "character",
    help = "Comma-delimited list of RDS objects containing count matrices"
  ),
  make_option(
    opt_str = "--output_name",
    type = "character",
    help = "Name to tag output matrix file with"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))
files <- strsplit(opts$matrix_files, ",")[[1]]

objlist = vector()
namelist = vector()
for (f in files) 
{
    pathparts <- strsplit(f, '/')[[1]]
    filename = pathparts[length(pathparts)]
    samplename = strsplit(filename, '\\.')[[1]][1]
    namelist <- append(namelist, samplename)
    countmatrix <- readRDS(f)
    seuratobj <- CreateSeuratObject(counts=countmatrix, project=samplename)
    objlist <- append(objlist, seuratobj)
}

if (length(objlist) == 1) {
    merged_obj <- objlist[[1]]
    } else if (length(objlist) == 2) {
    merged_obj <- merged_matrix <- merge(objlist[[1]], y=objlist[[-1]], add.cell.ids=namelist, project=opts$output_name)
    } else if (length(objlist) > 2) {
    merged_obj <- merge(objlist[[1]], y=objlist[-1], add.cell.ids=namelist, project=opts$output_name)
}

merged_matrix <- merged_obj[[]]
output_file <- paste(opts$output_name, '.merged.RDS', sep="")

saveRDS(merged_matrix, output_file)
