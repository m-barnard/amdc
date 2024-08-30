library(dplyr)
library(readr)
library(stringr)
library(rcompanion)
library(combinat)
library(clValid)
library(mclust)
library(parallel)
source('src/run_functions/new_mat_funcs.R')
source('src/run_functions/method_functions.R')
source('src/run_functions/clust_accuracy_functions.R')


run_sim <- function(file, dist_path, min_num, max_num,weight = FALSE, weight_val = NULL){
  df <- read_csv(file, show_col_types = FALSE)
  
  ### MY METHOD ### 
  if(weight == FALSE){
    indices <- list(c(1))
    type <- c('keep')
    trans_res <- final_gen(df$seqs, indices, type)
  } else if(weight == '9_5'){
    #pull 9-5 trans mat
    indices <- list(c(108, 204))
    type <- c('keep')
    trans_res_sub <- final_gen(df$seqs, indices, type)
    
    #pull full trans mat
    indices2 <- list(c(1))
    trans_res_full <- final_gen(df$seqs, indices2, type)
    
    #weight 9-5 
    w_vec <- trans_res_sub$vecs[[1]]*weight_val
    #get the non 9-5 part of the matrix
    sub_vec <- trans_res_full$vecs[[1]] - trans_res_sub$vecs[[1]]
    
    #get new set of trans vecs
    vec <- w_vec + sub_vec
    vec_s <- (vec*287)/unique(colSums(vec))
    
    #get new set of trans mats
    mat_s <- list()
    for(i in seq(1:ncol(vec_s))){
      mat_s[[i]] <- matrix(vec_s[,i],4,4)
    }
    trans_res <- list('vecs' = list('1_keep' = vec_s), 'mats' = list('1_keep' = mat_s))
  } else if (weight == 'O'){
    indices <- list(c(1))
    type <- c('keep')
    trans_res <- final_gen(df$seqs, indices, type)
    trans_res$vecs[[1]][6,] <- weight_val*trans_res$vecs[[1]][6,]
    trans_res$vecs[[1]] <- apply(trans_res$vecs[[1]], 2, function(x){287*x/sum(x)})
    get_mats <- list()
    for(i in seq(1, ncol(trans_res$vecs[[1]]))){
      get_mats[[i]] <- matrix(trans_res$vecs[[1]][,i],4,4)
    }
    trans_res$mats[[1]] <- get_mats
  }
  
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
  dist_mat <- read_rds(dist_path)
  
  

  df_list <- list()
  metric_df_list <- list()
  #lapply for running the k-means and hierarchical 10 times
  df_list <- mclapply(seq(1, 10), function(j){
    df_mat <- seq_cluster_options(df, ev_list,min_num,max_num, seed = j)
    df_clust <- seq_hcluster_options(df_mat, dist_mat, min_num, max_num, seed =j) %>%
      mutate(seed = j)
    return(df_clust)
  }, mc.cores = min(detectCores()[1], 10))
  metric_df_list <- mclapply(seq(1, 10), function(j){
    df_clust <- df_list[[j]]
    cluster_names <- colnames(df_clust)[!(colnames(df_clust) %in% c(colnames(df), 'seed'))]
    final_metric_df <- do.call(rbind,lapply(cluster_names, get_metrics_mat, df=df_clust, trans_res=trans_res, dist_mat = dist_mat))
    final_metric_df$seed <- j
    return(final_metric_df)
  }, mc.cores = min(detectCores()[1], 10))
  
  #combine the results from the lapply
  tot_metric <- do.call(rbind,metric_df_list)
  best_mat <- tot_metric %>%
    filter(type == 'mat') %>%
    group_by(cluster_num, type) %>%
    filter(metric == max(metric)) %>%
    ungroup()
  best_hier <- tot_metric %>%
    filter(type == 'hier') %>%
    group_by(cluster_num, type) %>%
    filter(dunn == max(dunn)) %>%
    ungroup()
  all_metric <- rbind(best_mat, best_hier)
    
    
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
      set.seed(5)
      best_seed <- sample(c(unique(sub$seed)), 1)
    }
    sub_best <- sub %>%
      filter(seed == best_seed) %>%
      dplyr::select(-seed)
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
  
  #save the final metric dataset
  final_metric_df <- do.call(rbind, final_metric_list) %>%
    arrange(type, cluster)
  
  return(list(final_metric_df, df_clust))
}