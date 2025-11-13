library(gemma.R)
library(argparse)

parser = ArgumentParser(description='Download DEA results from Gemma and save as TSV files.')
parser$add_argument('--username', type='character', help='Gemma username', default="raschwar")
parser$add_argument('--password', type='character', help='Gemma password', default="7nddtt")

parser$add_argument('--experiment', type='character',help='Gemma experiment accession', default="GSE280569")
args = parser$parse_args()
experiment <- args$experiment

dea_meta <- list()

set_gemma_user(args$username, args$password)

dea_meta <- get_dataset_differential_expression_analyses(experiment)
contrast_levels <- lapply(dea_meta[["experimental.factors"]], function(x) {
return(x[["value"]])
})
browser()
dea_meta[["level"]] <- sapply(contrast_levels, function(x) {
paste(x, collapse = ";")
})
# remove any columns with object data types
valid_columns <- sapply(dea_meta, function(col) {
!is.list(col) && !is.data.frame(col)
})

dea_meta <- dea_meta[, ..valid_columns]
# write to file
output_file <- paste0(experiment, "_meta.tsv")
write.table(dea_meta, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)