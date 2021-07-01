# R script to SoupX for removeal of mRNA contamination
# from single-cell data

options(error = traceback)

#install packages if not installed.
list.of.packages <- c("optparse", "codetools", "survival", "plyr", "zoo", "data.table", "htmlwidgets", "lazyeval", "reshape2", "SoupX")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
print(new.packages)
if(length(new.packages)) suppressMessages(install.packages(new.packages, repos='http://cran.us.r-project.org'))

suppressMessages(library(optparse))
suppressMessages(library(SoupX))

#proccess inputs
option_list <- list(
  make_option(
    opt_str = "--countmatrixdir",
    type = "character",
    help = "Path to directory containing raw and filtered matrix outputs from Cell Ranger count"
  ),
  make_option(
    opt_str = "--sample_name",
    type = "character",
    help = "Name associated with this sample"
  )
)

#parse options
opts <- parse_args(OptionParser(option_list = option_list))
countmatrixdir <- opts$countmatrixdir
sample_name <- opts$sample_name

# load 10X data
sc <- load10X(countmatrixdir)

# estimate contamination fraction rho
sc = autoEstCont(sc)

# clean the data
out_matrix = adjustCounts(sc)

output_file <- paste(sample_name, '.decontam.RDS', sep="")
saveRDS(out_matrix, output_file)
