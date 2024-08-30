library(dplyr)
library(stringr)
library(readr)
library(parallel)
source('src/gen_seqs_funcs/sim_seq_functions.R')
source('src/gen_seqs_funcs/get_TPMs.R')

#############################################
# Creates all simulated data for a given order of Markov chain
# i: sets the order of Markov chain to simulate from
# output: all data is saved as csv within simulated_sequences/order_{i}/
############################################

#keep this setting if running code on a slurm cluster, otherwise set i manually (1, 2, or 5 as options)
i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

P_all <- get_sim_TPMs(i)
start_state_list <- list('state_low' = list(1, 4), 'state_med' = list(1, 1, 4), 'state_high' = list(1,4, 1, 4),
                         'dur1_low' = list(4,4,4), 'dur1_med' = list(4,4,4), 'dur1_high' = list(4,4,4),
                         'dur2_low' = list(4,4,4), 'dur2_med' = list(4,4,4), 'dur2_high' = list(4,4,4),
                         'dur3_low' = list(4,4,4), 'dur3_med' = list(4,4,4), 'dur3_high' = list(4,4,4))
clust_prob_list <- list('state_low' = c(0.5, 0.5), 'state_med' = c(1/3, 1/3, 1/3), 'state_high' = rep(0.25, 4),
                        'dur1_low' = c(1/3, 1/3, 1/3), 'dur1_med' = c(1/3, 1/3, 1/3), 'dur1_high' = c(1/3, 1/3, 1/3),
                        'dur2_low' = c(1/3, 1/3, 1/3), 'dur2_med' = c(1/3, 1/3, 1/3), 'dur2_high' = c(1/3, 1/3, 1/3),
                        'dur3_low' = c(1/3, 1/3, 1/3), 'dur3_med' = c(1/3, 1/3, 1/3), 'dur3_high' = c(1/3, 1/3, 1/3))

folders <- names(start_state_list)
for(f in folders){
  print(f)
  #create directory to save the results
  dir.create(file.path(paste0('simulated_sequences/order_', i), f))
  #seq(0,499) defines the # of datasets, here we are creating 500
  null_res <- mclapply(seq(0, 499), function(k){
    print(k)
    #250 in the seq() function below defines the # of sequences per dataset
    seq_data <- as.data.frame(t(sapply(seq(1 + 250*k,250 + 250*k), gen_clust_seq, P_list = P_all[[f]], order = i, 
                                       num.iters = 500, start_state = start_state_list[[f]], clust_prob = clust_prob_list[[f]])))
   #write the simulated dataset to the folder we created above
    write_csv(seq_data, paste0('simulated_sequences/order_', i, '/', f, '/', k, '.csv'))
    invisible()
  }, mc.cores = 10)
}

