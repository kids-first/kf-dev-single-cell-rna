# R script for mRNA decontamination of 10X single-cell data using SoupX

options(error = traceback)

#install packages if not installed.
list.of.packages <- c("optparse", "codetools", "survival", "plyr", "zoo", "data.table", "htmlwidgets", "lazyeval", "reshape2", "SoupX")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(new.packages)
if(length(new.packages)) suppressMessages(install.packages(new.packages, repos='http://cran.us.r-project.org'))

suppressMessages(library(optparse))
suppressMessages(library(SoupX))
suppressMessages(library(Seurat))
suppressMessages(library(DropletUtils))

#process inputs
option_list <- list(
  make_option(
    opt_str = "--raw",
    type = "character",
    help = "Path to h5 file containing raw matrix from Cell Ranger count"
  ),
  make_option(
    opt_str = "--fil",
    type = "character",
    help = "Path to h5 file containing filtered matrix from Cell Ranger count"
  ),
  make_option(
    opt_str = "--cluster",
    type = "character",
    help = "Path to csv file containing barcodes and cluster id from Cell Ranger count"
  ),
  make_option(
    opt_str = "--sample_name",
    type = "character",
    help = "Name associated with this sample"
  )
)

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
raw <- opts$raw
fil <- opts$fil
clusters_file <- opts$cluster
sample_name <- opts$sample_name

# load 10X data and convert to Soup Channel (object SoupX uses for analysis)
fil_mat <- Read10X_h5(fil, use.names = T)
raw_mat <- Read10X_h5(raw, use.names = T)
clusters <- read.csv(clusters_file)
mDat <- data.frame(clusters=clusters$Cluster,row.names=clusters$Barcode)
sc <- SoupChannel(raw_mat, fil_mat, mDat)

# estimate contamination fraction rho
sc <- autoEstCont(sc)

# clean the data
out_matrix <- adjustCounts(sc, roundToInt = T)
colnames(out_matrix) = gsub('-1$', '', colnames(out_matrix))
colnames(out_matrix) = paste(sample_name, colnames(out_matrix), sep=":")

saveRDS(out_matrix, paste0(sample_name, ".rds"))

DropletUtils:::write10xCounts(sample_name, out_matrix)