suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))
source('src/run_functions/calc_sim_clust.R')


######################
# Run our method and hierarchical clustering and assess cluster accuracy for all datasets in a simulation scenario
# NOTE: simulation inputs are mainly created through reading in folder names for the generated simulated sequences
# so the inputs may require some tweaking if all sequences are not simulated before this
# i: indexes the input_grid matrix to identify the simulation scenario to run
# output: all data is saved as rds within results/order_{}/ and named by simulation scenario
######################

#keep this setting if running code on a slurm cluster, otherwise set i manually 1-36
i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#create inputs 
outer_folders <- sort(list.files('simulated_sequences/'))
folders <- sort(list.files('simulated_sequences/order_1/'))

input_grid <- expand.grid('o_fold' = outer_folders, 'folder_name'  = folders)
input_grid <- cbind(input_grid, 'max' = ifelse(input_grid[,'folder_name' ] == 'state_high', 7, 
                                       ifelse(input_grid[,'folder_name' ] == 'state_low', 5, 6)))
inputs <- input_grid[i,]
print(inputs)

input_path <- paste0('simulated_sequences/',inputs$o_fold,'/', inputs$folder_name, '/')

#run simulation on each of the 500 datasets - change seeds assignment if you want to use fewer/more datasets
seeds <- seq(0, 499)
out <- mclapply(seeds, run_sim, path = input_path, min_num = 2, max_num = inputs$max, mc.cores = 24)

clust_acc <- do.call(rbind, lapply(out, function(x){x[[1]]}))
clust_metric <- do.call(rbind, lapply(out, function(x){x[[2]]}))


#save results for all datasets within a given simulation scenario
write_rds(clust_acc, paste0('results/', inputs$o_fold,'/', inputs$folder_name, '_cluster_accuracy.rds'))
write_rds(clust_metric, paste0('results/',inputs$o_fold,'/', inputs$folder_name, '_cluster_metrics.rds'))
print(paste0('Done with ', inputs$o_fold,'/', inputs$folder_name))
