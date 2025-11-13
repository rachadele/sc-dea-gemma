library(argparse)
library(tidyr)
parser <- ArgumentParser(description='Convert RDS DEA results to TSV format')
parser$add_argument('--dea_meta', type='character', help='Path to the input RDS meta file', default = "/space/scratch/DEA_testing/class_results/class_meta.rds")

args <- parser$parse_args()
# read in RDS files
dea_meta <- readRDS(args$dea_meta)

experiments_meta <- names(dea_meta)

experiments <- names(dea_meta)
# initialize list for each experiment

for (experiment in experiments) {
  message("Processing experiment: ", experiment)
  
  meta_table <- dea_meta[[experiment]]
  # extract values associated with each experimental factor
  contrast_levels <- lapply(meta_table[["experimental.factors"]], function(x) {
    return(x[["value"]])
  })
  browser()
  meta_table[["level"]] <- sapply(contrast_levels, function(x) {
    paste(x, collapse = ";")
  })
  # remove any columns with object data types
  valid_columns <- sapply(meta_table, function(col) {
    !is.list(col) && !is.data.frame(col)
  })

  meta_table <- meta_table[, ..valid_columns]
  # write to file
  output_file <- paste0(experiment, "_meta.tsv")
  write.table(meta_table, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
  


