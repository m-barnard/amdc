library(readr)
library(dplyr)
library(parallel)
source('src/cluster_stability/run_one_clust.R')

get_boot_clust <- function(seed, df, df_full,comp_clust_df, opt = T){
  #############################
  # Get jaccard and clustering result for bootstrapped dataset
  # df: original data
  # df_full: original data + eigenvectors from SVD
  # min_num: minimum cluster number to run
  # max_num: maximum cluster number to run
  # opt: T indicates we will select the # of clusters with our metric (this is what the main results are from)
  #      F indicates that we select 8 clusters every time as that was true of the originla data
  # outputs: dataset with the jaccard index for each observation, best # of clusters for the bootstrapped dataset
  ############################
  set.seed(seed)
  indvs <- sample(unique(df$user_id), length(unique(df$user_id)), replace = T)
  boot_df <- do.call(rbind, lapply(seq(1, length(indvs)), function(i){
    set.seed(seed + i)
    k <- indvs[i]
    sub_df <- df %>%
      filter(user_id == k) 
    sub_df$n <- nrow(sub_df)
    sub_boot_df <- sub_df %>%
      sample_n(size = unique(n), replace = T)
    return(sub_boot_df)
  }))
  
  ## get clusters using bootstrapped dataset ##
  out_boot <- run_sim_mat(boot_df, 2, 9)
  
  ### select best cluster by our metric or just select 8 ###
  if(opt == T){
    best_clust_boot <- out_boot[[1]] %>%
      filter(metric == max(metric)) %>%
      pull(cluster)
    best_seed_boot <- out_boot[[1]]%>%
      filter(metric == max(metric)) %>%
      pull(seed)
  } else if(opt == F){
    best_clust_boot <- colnames(out_boot[[2]])[str_detect(colnames(out_boot[[2]]), 'clust8')]
    best_seed_boot <- out_boot[[1]] %>%
      filter(cluster == best_clust_boot) %>%
      pull(seed)
  }
  ### get the new clusters based on the centers of the bootstrapped data clusters ###
  full_boot_df <- out_boot[[3]][[best_seed_boot]]
  best_center_cols <- colnames(full_boot_df)[str_detect(colnames(full_boot_df), best_clust_boot)]
  rename_cols <- ifelse(str_detect(best_center_cols, 'center'), str_sub(best_center_cols, nchar(best_center_cols), nchar(best_center_cols)), best_center_cols)
  center_df <- full_boot_df %>%
    dplyr::select(all_of(best_center_cols)) %>%
    arrange(!!sym(best_clust_boot)) %>%
    unique()
  colnames(center_df) <- rename_cols
  num_evs <- ncol(center_df) - 1
  min_dist_boot_clust <- unlist(mclapply(seq(1, nrow(df_full)), function(i){
    evs <- df_full[i, as.character(seq(1, num_evs))]
    dist_vec <- sapply(seq(1, nrow(center_df)), function(j){
      return(norm_vec(unlist(center_df[j, as.character(seq(1, num_evs))] - evs)))
    })
    boot_clust <- which(dist_vec == min(dist_vec))
    return(boot_clust)
    }, mc.cores = 5))
  clust_df <- data.frame('uniq_id' = df_full$uniq_id, 'cluster' = min_dist_boot_clust)
  
  ### calc jaccard index ###
  jac_vec <- unlist(mclapply(seq(1, length(clust_df$uniq_id)), function(x){
    boot_x <- clust_df %>%
      filter(cluster == clust_df[x, 'cluster'])
    comp_x_clust <- comp_clust_df %>%
      filter(uniq_id == clust_df$uniq_id[x]) %>%
      pull(clust8_EV1234)
    comp_x <- comp_clust_df %>%
      filter(clust8_EV1234 == comp_x_clust)
    return(jaccard(boot_x$uniq_id, comp_x$uniq_id))
  }, mc.cores = 5))
  final_df <- data.frame('uniq_id' = clust_df$uniq_id, 'jac' = jac_vec, 'best_clust' = best_clust_boot, 'seed' = seed)
  return(final_df)
}
