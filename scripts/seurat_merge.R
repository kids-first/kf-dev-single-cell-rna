# R script for merging 10X single-cell datasets using Seurat

options(error = traceback)

suppressMessages(library(optparse))
suppressMessages(library(Seurat))

#process inputs
option_list <- list(
  make_option(
    opt_str = "--matrix_dirs",
    type = "character",
    help = "Comma-delimited list of RDS objects containing count matrices"
  ),
  make_option(
    opt_str = "--doublets_files",
    type = "character",
    help = "Csv files with doublet information"
  ),
  make_option(
    opt_str = "--align_qc_files",
    type = "character",
    help = "QC RDS file from D3b alignment workflow"
  ),
  make_option(
    opt_str = "--output_name",
    type = "character",
    help = "Name to tag output matrix file with"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))
mats <- strsplit(opts$matrix_dirs, ",")[[1]]
doubs <- strsplit(opts$doublets_files, ",")[[1]]
align_qcs <- strsplit(opts$align_qc_files, ",")[[1]]

objlist = vector()
namelist = vector()
i <- 1
for (m in mats)
{
    #figure out sample name from file basename
    pathparts <- strsplit(m, '/')[[1]]
    filename = pathparts[length(pathparts)]
    samplename = strsplit(filename, '\\.')[[1]][1]

    #add name to list of names
    namelist <- append(namelist, samplename)

    #read in matrix and make seurat object
    countmatrix <- Read10X(m)
    seuratobj <- CreateSeuratObject(counts=countmatrix, project=samplename)

    #read doublet file and remove doublets
    doublet_file <- doubs[[i]]
    doublets <- read.table(doublet_file, sep = ",", , header=F, row.names=1)
    colnames(doublets) <- c("Doublet_score","Is_doublet")
    seuratobj <- AddMetaData(seuratobj,doublets)

    #mark barcodes where Is_doublet column is False as QC Passes
    seuratobj[['QC']] <- ifelse(
    seuratobj@meta.data$Is_doublet == 'True','Doublet','Pass')

    # read align QC file and merge
    align_qc_fn <- align_qcs[[i]]
    align_qc = readRDS(align_qc_fn)
    seuratobj <- merge(seuratobj, align_qc)

    #add subset of QC passes list of objects
    objlist <- append(objlist, subset(seuratobj, subset = QC == 'Pass'))
    i <- i + 1
}

if (length(objlist) == 1) {
    merged_obj <- objlist[[1]]
    } else if (length(objlist) == 2) {
    merged_obj <- merged_matrix <- merge(objlist[[1]], y=objlist[[-1]], add.cell.ids=namelist, project=opts$output_name)
    } else if (length(objlist) > 2) {
    merged_obj <- merge(objlist[[1]], y=objlist[-1], add.cell.ids=namelist, project=opts$output_name)
}

#export merged matrix and Seurat object
merged_matrix <- merged_obj[[]]
output_file <- paste(opts$output_name, '.merged_matrix.RDS', sep="")
saveRDS(merged_matrix, output_file)

output_file <- paste(opts$output_name, '.merged_object.RDS', sep="")
saveRDS(merged_obj, output_file)
