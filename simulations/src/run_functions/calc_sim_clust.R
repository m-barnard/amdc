suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rcompanion))
suppressPackageStartupMessages(library(combinat))
suppressPackageStartupMessages(library(clValid))
suppressPackageStartupMessages(library(mclust))
source('src/run_functions/new_mat_funcs.R')
source('src/run_functions/method_functions.R')
source('src/run_functions/clust_accuracy_functions.R')


run_sim <- function(seed, path, min_num, max_num){
  #############################
  # Run our method and hierarchical clustering and assess cluster accuracy for a single simulated dataset
  # path: path to the simulated data folder
  # min_num: minimum cluster number to run
  # max_num: maximum cluster number to run
  # outputs: list of datasets; first contains cluster accuracy results, 
  #                            second contains the cluster evaluation metrics for all cluster results
  ############################
  set.seed(seed)
  print(seed)
  df <- read_csv(paste0(path, seed, '.csv'), show_col_types = FALSE)
  num_clusts <- length(unique(df$cluster))
  
  ### Prep data for our method ### 
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
  df <- cbind(df, svd_res$v[,seq(1,max_evs)])
  colnames(df) <- c(og_names, as.character(seq(1, max_evs)))
  ev_list <- lapply(seq(1,max_evs), function(x){seq(1,x)})
  
  ### Get distance matrix for hierarhcical clustering ### 
  dist_mat <- adist(df$seqs)
  
  
  ### Run kmeans and hierarhcical clustering 10 times for each cluster # and # of EVs### 
  df_list <- list()
  metric_df_list <- list()
  for(k in 1:10){
    df_mat <- seq_cluster_options(df, ev_list,min_num,max_num, seed = k)
    df_clust <- seq_hcluster_options(df_mat, dist_mat, min_num, max_num, seed =k)
    cluster_names <- colnames(df_clust)[!(colnames(df_clust) %in% colnames(df))]
    final_metric_df <- do.call(rbind,lapply(cluster_names, get_metrics_mat, df=df_clust, trans_res=trans_res, dist_mat = dist_mat))
    final_metric_df$seed <- k
    df_clust$seed <- k
    df_list[[k]] <- df_clust
    metric_df_list[[k]] <- final_metric_df
  }
  #combine the results from the loops
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
      filter((cluster_num == index[k, 'cluster_num']) & (type == index[k, 'type']))
      #to get rid of weird sample behavior when x has length 1
      if(length(unique(sub$seed)) ==1){
        best_seed <- unique(sub$seed)[1]
      } else{
        best_seed <- sample(unique(sub$seed), 1)
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
  #save the final metric dataset, df_clust is already saved
  final_metric_df <- do.call(rbind, final_metric_list) %>%
    arrange(type, cluster)
  
  mat_col <- final_metric_df %>%
    filter(cluster_num == num_clusts & type == 'mat') %>%
    pull(cluster)
  hier_col <- final_metric_df %>%
    filter(cluster_num == num_clusts & type == 'hier') %>%
    pull(cluster)
  
  
  best_mat <- final_metric_df %>%
    filter(type == 'mat') %>%
    filter(metric == max(metric))
  #fix if two clusterings have the same metric (unlikely but could happen)
  if(nrow(best_mat) > 1){
    best_mat <- best_mat %>%
      filter(cluster_num == min(cluster_num))
  }
  
  best_hier <- final_metric_df %>%
    filter(type == 'hier') %>%
    filter(dunn == max(dunn))
  #fix if two clusterings have the same metric (unlikely but could happen)
  if(nrow(best_hier) > 1){
    best_hier <- best_hier %>%
      filter(cluster_num == min(cluster_num))
  }
  
  ### Calculate cluster evaluation metrics ###
  #contingency tables
  mat_table <- table(as.factor(df_clust$cluster), as.factor(df_clust[,mat_col]))
  h_table <- table(as.factor(df_clust$cluster), as.factor(df_clust[,hier_col]))

  #get the ARI for the correct # of clusters
  ARI_h <- adjustedRandIndex(as.factor(df_clust$cluster), as.factor(df_clust[,hier_col]))
  ARI_m <- adjustedRandIndex(as.factor(df_clust$cluster), as.factor(df_clust[,mat_col]))
  
  #get the ARI for the best # of clusters
  ARI_hbest <- adjustedRandIndex(as.factor(df_clust$cluster), as.factor(df_clust[,best_hier$cluster]))
  ARI_mbest <- adjustedRandIndex(as.factor(df_clust$cluster), as.factor(df_clust[,best_mat$cluster]))
  
  out <- data.frame(rbind(c(calc_accuracy(mat_table, num_clusts), 'ARI' = ARI_m, 'ARI_best' = ARI_mbest, 'best_clust' = best_mat$cluster_num), 
                          c(calc_accuracy(h_table, num_clusts), 'ARI' = ARI_h, 'ARI_best' = ARI_hbest, 'best_clust' = best_hier$cluster_num))) %>%
    mutate(type = c('mat', 'hier')) %>%
    mutate(seed = seed)
  final_metric_df <- final_metric_df %>%
    mutate(seed = seed)
  print(paste0('done ', seed))
  return(list(out, final_metric_df))
}
  













