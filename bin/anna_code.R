library(gemma.R)

sc_experiment_list <- c('GSE267301',
                        'GSE213364',
                        'GSE157827',
                        'GSE163323.3',
                        'GSE133283',
                        'GSE221712',
                        'GSE183757',
                        'GSE283187',
                        'GSE181786',
                        'GSE280569')
dea_res <- list()
for (experiment in sc_experiment_list){
  sc_DE <- get_differential_expression_values(experiment)
  dea_res[[experiment]] <- sc_DE
}

dea_meta <- list()
for (experiment in sc_experiment_list){
  sc_meta <- get_dataset_differential_expression_analyses(experiment)
  dea_meta[[experiment]] <- sc_meta
}