suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(combinat))
suppressPackageStartupMessages(library(parallel))
source('src/run_functions/clust_accuracy_functions.R')


######################
# Get accuracy of nTreeClus clusters
# NOTE: simulation inputs are mainly created through reading in folder names for the generated simulated sequences
# so the inputs may require some tweaking if all sequences are not simulated before this
# output: all data is saved as rds within results/order_{}/ and named by simulation scenario
######################

outer_folders <- sort(list.files('simulated_sequences/'))
folders <- sort(list.files('simulated_sequences/order_1/'))
input_grid <- expand.grid('o_fold' = outer_folders, 'folder_name'  = folders)

#i indexes the input grid
for(i in seq(1, 36)){
  inputs <- input_grid[i,]
  res_df <- read_csv(paste0('results/', inputs$o_fold, '/', inputs$folder_name, '_ntrees_clust.csv'))
  
  #0 - 499 indexes the simulated datasets - change if you want fewer/more datasets
  df_list <- lapply(seq(0, 499), function(i){
    return(res_df %>% filter(seed == i))
  })
  
  clean_res <- mclapply(df_list, get_mult_acc, mc.cores = 10)
  clean_res_df <- as.data.frame(do.call(rbind, clean_res))
  write_rds(clean_res_df, paste0('results/', inputs$o_fold, '/', inputs$folder_name, '_ntrees_accuracy.csv'))
}
