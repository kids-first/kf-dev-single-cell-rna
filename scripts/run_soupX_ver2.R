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

# get contamination summary table per library
contamination_summary <- matrix(ncol=4, nrow=1) %>% as.data.frame()
colnames(contamination_summary) <- c("library", "SoupX_contamination", "adjusted", "adjusted_SoupX_cont")

# load 10X data and convert to Soup Channel (object SoupX uses for analysis)
fil_mat <- Read10X_h5(fil, use.names = T)
raw_mat <- Read10X_h5(raw, use.names = T)
clusters <- read.csv(clusters_file)
mDat <- data.frame(clusters=clusters$Cluster,row.names=clusters$Barcode)
sc <- SoupChannel(raw_mat, fil_mat, mDat)

# estimate contamination fraction rho
# It can be higheer in some cases, but we shouldn't remove more than 20% background
sc1 = try(autoEstCont(sc))

    if (class(sc1)=="try-error") {
        cat("Estimated soup fraction >50%, set to 5% by default.\n")
        sc1 = setContaminationFraction(sc, 0.05)
        contamination_summary[i,] <- c(sample_name, 0.5, "yes", 0.05)
    } else if (sc1$fit$rhoEst>0.2) {
        print(sc1$fit$rhoEst)
        cat("Estimated soup fraction >20%, set to 5% by default.\n")
        contamination_summary[i,] <- c(sample_name, sc1$fit$rhoEst, "yes", 0.05)
        sc1 = setContaminationFraction(sc, 0.05)
    } else {
        contamination_summary[i,] <- c(sample_name, sc1$fit$rhoEst, "no", sc1$fit$rhoEst)
    }

#Clean the data
out = adjustCounts(sc1, roundToInt = T)
colnames(out) <- paste0(sample_name, ":", colnames(out))

# save objects
saveRDS(sc1, paste0(output_dir, "sc_soupX_", sample_name,".rds"))
saveRDS(out, paste0(output_dir, "soupX_count_mtx_", sample_name,".rds"))

print("Merge mtxs")
  if(i==1){
        out_merge <- out
    } else{
        out_merge <- cbind(out_merge, out)
    }


DropletUtils:::write10xCounts(output_dir, out_merge)


# save combined filtered matrix and 
saveRDS(out_merge, paste0(out_path, "soupX_count_mtx_combined.rds"))

# save contamination summary table
write.csv(contamination_summary, paste0(out_path,"soupX_contamination_summary.csv"))
