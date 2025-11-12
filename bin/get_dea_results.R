library(gemma.R)
library(argparse)

parser = ArgumentParser(description='Download DEA results from Gemma and save as TSV files.')
parser$add_argument('--username', type='character', help='Gemma username', required=TRUE)
parser$add_argument('--password', type='character', help='Gemma password', required=TRUE)

parser$add_argument('--experiment', type='character',
					help='Gemma experiment accession', required=TRUE)
args = parser$parse_args()
experiment <- args$experiment

dea_res <- list()

set_gemma_user(args$username, args$password)

sc_DE <- get_differential_expression_values(experiment, readableContrasts = FALSE)
for (resultSet in names(sc_DE)) {
	# save tsv with resultSet and experiment in filename
	filename <- paste0(experiment, "_", resultSet, "_dea_results.tsv")
	write.table(sc_DE[[resultSet]], file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}