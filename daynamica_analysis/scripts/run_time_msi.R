suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rcompanion))
suppressPackageStartupMessages(library(combinat))
suppressPackageStartupMessages(library(clValid))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rbenchmark))
source('src/run_functions/new_mat_funcs.R')
source('src/run_functions/method_functions.R')
source('src/run_functions/clust_accuracy_functions.R')
source('src/run_functions/run_clust_time.R')


### HIER TIME ###
out0 <- benchmark(run_df_prep_hier('clean_data/clean_day_seqs.csv'), replications = 10)
write_rds(out0, 'results_time/prep_day_hier.rds')

out02 <- benchmark(run_df_prep_hier('clean_data/clean_week_seqs.csv'), replications = 50)
write_rds(out02, 'results_time/prep_week_hier.rds')

dist_mat <- read_rds('clean_data/clean_day_distmat.rds')
out <- benchmark(run_sim_hier('clean_data/clean_day_seqs.csv', 2, 9, dist_mat))
write_rds(out, 'results_time/clust_day_hier.rds')

dist_mat <- read_rds('clean_data/clean_week_distmat.rds')
out2 <- benchmark(run_sim_hier('clean_data/clean_week_seqs.csv', 2, 9, dist_mat))
write_rds(out2, 'results_time/clust_week_hier.rds')

### MAT TIME ###
mout0 <- benchmark(run_df_prep_mat('clean_data/clean_day_seqs.csv'))
write_rds(mout0, 'results_time/prep_day_mat.rds')

mout02 <- benchmark(run_df_prep_mat('clean_data/clean_week_seqs.csv'))
write_rds(mout02, 'results_time/prep_week_mat.rds')

mout <- benchmark(run_sim_mat('clean_data/clean_day_seqs.csv', 2, 9))
write_rds(mout, 'results_time/clust_day_mat.rds')

mout2 <- benchmark(run_sim_mat('clean_data/clean_week_seqs.csv', 2, 9))
write_rds(mout2, 'results_time/clust_week_mat.rds')







