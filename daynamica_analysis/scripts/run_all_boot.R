library(readr)
library(dplyr)
library(tibble)
library(parallel)
source('src/cluster_stability/run_boot_proc.R')

input_grid <- expand.grid(seq(1,500), c(T))
names(input_grid) <- c('seed', 'opt')

val <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #input number 1-20 if not using the bash script
all_nums <- seq(1 + (val-1)*25, val*25)

df <- read_csv('clean_data/clean_day_seqs.csv', show_col_types = FALSE)

#get full df with the EVs
indices <- list(c(1))
type <- c('keep')
trans_res <- final_gen(df$seqs, indices, type)
#run first svd step
trans_mat <- trans_res[["vecs"]]$`1`
trans_mat_s <- scale(trans_mat, scale = FALSE)
svd_res <- svd(trans_mat_s)
perc_var <- sapply(seq(1, length(svd_res$d)), calc_perc, vec = svd_res$d, n = 16)
#only using eigenvectors explaining greater than 1% of variance
max_evs <- length(perc_var[perc_var > 0.01])
og_names <- colnames(df)
df_full <- cbind(df, svd_res$v[,seq(1,max_evs)])
colnames(df_full) <- c(og_names, as.character(seq(1, max_evs)))

#read in cluster results on the original data
comp_df <- read_csv('results/amdc_clust_res.csv', show_col_types = FALSE)

#run the stability procedure for 25 bootstrapped datasets
#NOTE: if you do not want to use multi-node parallelization here, instead of all_nums just use seq(1,500)
out <- mclapply(all_nums, function(x){
  print(x)
  final_df <- get_boot_clust(input_grid$seed[x], df = df, df_full = df_full, comp_clust_df = comp_df, opt = input_grid$opt[x])
}, mc.cores = 5)
out_df <- do.call(rbind, out)

write_rds(out_df, paste0('results_stability/', val, '.rds'))
