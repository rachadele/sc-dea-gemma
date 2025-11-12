library(argparse)
library(tidyr)
library(dplyr)
parser <- ArgumentParser(description='Convert RDS DEA results to TSV format')
parser$add_argument('--dea_results', type='character', help='Path to the input RDS results file', default = "/space/scratch/DEA_testing/class_results/class_res.rds")


args <- parser$parse_args()
dea_results <- readRDS(args$dea_results)

experiments <- names(dea_results)

for (exp in experiments) {
  res <- dea_results[[exp]]

  contrasts <- names(res)
  for (ct in contrasts) {
	ct_df <- as.data.frame(res[[ct]])
	if (nrow(ct_df) == 0) 
	next
	ct_df$contrast <- ct
	browser()
	output_file <- paste0(exp, "_", ct, "_de_res.tsv")
	write.table(ct_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}