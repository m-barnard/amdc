library(readr)
library(dplyr)
library(stringr)
library(tibble)
source('src/run_functions/new_mat_funcs.R')
source('src/run_functions/method_functions.R')
source('src/run_functions/clust_accuracy_functions.R')


run_sim_mat <- function(df, min_num, max_num){
  #############################
  # Run our method and and assess cluster accuracy for a single bootrapped dataset
  # df: bootstrapped dataset
  # min_num: minimum cluster number to run
  # max_num: maximum cluster number to run
  # outputs: list of datasets; first contains cluster accuracy results, 
  #                            second contains the cluster evaluation metrics for all cluster results
  #                            third contains all cluster results and the centers for all clusters
  ############################
  indices <- list(c(1))
  type <- c('keep')
  trans_res <- final_gen(df$seqs, indices, type)
  #run first svd step
  trans_mat <- trans_res[["vecs"]]$`1`
  trans_mat_s <- scale(trans_mat, scale = FALSE)
  svd_res <- svd(trans_mat_s)
  
  
  
  #run output with exact # of clusters
  perc_var <- sapply(seq(1, length(svd_res$d)), calc_perc, vec = svd_res$d, n = 16)
  #only using eigenvectors explaining greater than 1% of variance
  max_evs <- length(perc_var[perc_var > 0.01])
  og_names <- colnames(df)
  df <- cbind(df, svd_res$v[,seq(1,max_evs)])
  colnames(df) <- c(og_names, as.character(seq(1, max_evs)))
  ev_list <- lapply(seq(1,max_evs), function(x){seq(1,x)})
  
  
  # run k means 10 times
  df_full_list <- mclapply(seq(1, 10), function(j){
    df_clust <- seq_cluster_options_centers(df, ev_list,min_num,max_num, seed = j) %>%
      mutate(seed = j)
    return(df_clust)
  }, mc.cores = 10)
  df_list <- lapply(df_full_list, function(x_df){
    select_cols <- colnames(x_df)[str_detect(colnames(x_df), 'center') == F]
    small_df <- x_df %>%
      dplyr::select(all_of(select_cols))
    return(small_df)
  })
  metric_df_list <- mclapply(seq(1, 10), function(j){
    df_clust <- df_list[[j]]
    cluster_names <- colnames(df_clust)[!(colnames(df_clust) %in% c(colnames(df), 'seed'))]
    final_metric_df <- do.call(rbind,lapply(cluster_names, get_metrics_mat, df=df_clust, trans_res=trans_res, dist_mat = dist_mat))
    final_metric_df$seed <- j
    return(final_metric_df)
  }, mc.cores = 10)
  
  tot_metric <- do.call(rbind,metric_df_list)
  all_metric <- tot_metric %>%
      filter(type == 'mat') %>%
      group_by(cluster_num, type) %>%
      filter(metric == max(metric)) %>%
      ungroup()
    
    
  all_df <- do.call(rbind, df_list)
    
  index <- all_metric %>%
    dplyr::select(cluster_num, type) %>%
    unique() %>%
    as.matrix()
    
  df_clust <- df
  final_metric_list <- list()
  #get the result for each cluster # and method (because some of the 10 runs give identical clustering results)
  for (k in seq(1, nrow(index))){
    sub <- all_metric %>%
      filter((cluster_num == as.vector(index[k, 'cluster_num'])) & (type == as.vector(index[k, 'type'])))
    #to get rid of weird sample behavior when x has length 1
    if(length(unique(sub$seed)) ==1){
      best_seed <- unique(sub$seed)[1]
    } else{
      best_seed <- sample(c(unique(sub$seed)), 1)
    }
    sub_best <- sub %>%
      filter(seed == best_seed)
    if(nrow(sub_best) > 1){
        sub_best <- sub_best %>%
          filter(evs == min(evs))
      }
    sub_cluster_name <- sub_best$cluster

    sub_df <- all_df %>%
      filter(seed == best_seed)
    df_clust <- cbind(df_clust, sub_df[sub_cluster_name])
    final_metric_list[[k]] <- sub_best
  }
  #save the final metric dataset, df_clust is already saved
  final_metric_df <- do.call(rbind, final_metric_list) %>%
    arrange(type, cluster)
    
    
  
  return(list(final_metric_df, df_clust, df_full_list))
}

jaccard <- function(a, b) {
  ####################
  # Calculate the jaccard index for two sets a,b
  ####################
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


norm_vec <- function(x){ sqrt(sum(x^2)) }
