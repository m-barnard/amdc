##### Helper functions to run the method ###
library(tibble)

### Get % of variance explained by the principal axes ###
calc_perc <- function(i, vec, n){
  ###################
  # vec: vector of singular values from SVD
  ##################
  #this is the relationship between the singular values and the eigenvalues
  return((vec[i]^2/(n-1))/sum(vec^2/(n-1)))
}

### Wrapper fucntions for stanndard clustering algorithms ###
kmeans_mb <- function(cols, df, k, seed){
  set.seed(seed)
  return(kmeans(df[cols], k)$cluster)
}
hclust_mb <- function(dist_mat, seed){
  set.seed(seed)
  return(hclust(as.dist(dist_mat), method = 'average'))
}

### Wrapper function to run kmeans given a set of eigenvectors ###
clust2 <- function(evs, clust_num, df, seed){
  ########
  # evs: vector of eigenvector numbers (e.g., c(1,2) for the first two eigenvectors)
  ########
  evs_char <- as.character(evs)
  col_name <- paste0('clust', clust_num, '_', 'EV', paste0(evs, collapse = ''))
  out_df <- data.frame(kmeans_mb(evs_char, df, clust_num, seed))
  colnames(out_df) <- c(col_name)
  return(out_df)
}
### Wrapper function to run hierarchical clustering for a certain cluster ####
hclust_func <- function(clust_num, dist_mat, seed){
  clust_obj <- hclust_mb(dist_mat, seed)
  out_df <- data.frame(cutree(clust_obj, clust_num))
  colnames(out_df) <- paste0('clust', clust_num)
  return(out_df)
}


### Cluster evaluation metrics functions ###
comp_norm <- function(mat_i, svd_avg){
  #######################
  # Calculate the metrics proposed in equations 1 and 2
  # mat_i: a single adjacency matrix
  # svd_avg: SVD results object for an average adjacency matrix
  #####################
  res <- svd(mat_i)
  approx_i <- svd_avg$u %*% diag(res$d) %*% t(svd_avg$v)
  diff <- mat_i - approx_i
  return(norm(diff, type  ='F')^2)
}

get_cluster_info <- function(cluster, df, trans_res, clust_name){
  ################
  # Calculate metrics (equation 1) for all sequences in a cluster
  # cluster: integer, clsuter #
  # trans_res: results object produced by final_gen() (set of ajdacency matrices)
  # clust_name: column name that contains the clustering output
  ################
  trans_all_mat <- trans_res[["mats"]]$`1`
  bool_vec = sapply(df[clust_name], function(x){x ==cluster})
  index_cols <- which(bool_vec == TRUE)
  clust_mats <- trans_all_mat[index_cols]
  avg_clust <- Reduce('+', clust_mats)/length(clust_mats)
  svd_a <- svd(avg_clust)
  clust_norms <- sum(sapply(clust_mats, comp_norm, svd_avg = svd_a))
  return(list('mat' = avg_clust, 'norm' = clust_norms))
}

get_metrics_mat <- function(cluster_name, df, trans_res, dist_mat = NULL){
  ######################
  # Get all cluster evaluation metrics for a given cluster output
  # cluster_name: column name that contains the clustering output
  # trans_res: results object produced by final_gen() (set of adjacency matrices)
  ######################
  clust_num <- as.numeric(str_sub(cluster_name,6,6)) ##NOTE: CAN'T DO MORE THAN 9 CLUSTERS WITH THIS
  if(str_detect(cluster_name, '_')){
    evs <- str_sub(cluster_name,-1)
    type <- 'mat'
    dunn <- NA
  } else{
    evs <- NA
    type <- 'hier'
    dunn <- dunn(dist_mat, df[cluster_name])
  }
  
  ### WITHIN CALC ###
  cluster_data <- lapply(seq(1, clust_num), get_cluster_info, df = df, trans_res = trans_res, clust_name = cluster_name)
  within <- sum(sapply(cluster_data, function(x){x$'norm'}))
  avg_mats <- lapply(cluster_data, function(x){x$'mat'})
  
  ### BETWEEN CALC ###
  trans_all_mat <- trans_res[["mats"]]$`1`
  o_avg_clust <- Reduce('+', trans_all_mat)/length(trans_all_mat)
  svd_o_avg <- svd(o_avg_clust)
  between_norms <- sapply(avg_mats, comp_norm, svd_avg = svd_o_avg)
  n <- df %>%
    group_by(!!sym(cluster_name)) %>%
    summarize(n = n()) %>%
    arrange(!!sym(cluster_name)) #arranges into same numerical order as avg_mats, and therefore between_norms
  n_vec <- n$n
  between <- sum(n_vec*between_norms)
  
  metric <- between*(sum(n_vec) - clust_num)/(within*(clust_num-1))
  
  return(data.frame('var_with' = within, 'var_between'= between, 'metric'= metric, 'dunn' = dunn, 'cluster' = cluster_name, 'cluster_num' = as.character(clust_num), 'evs'= evs, 'type'=type))
}

### Wrapper functions for getting a clustering results for a set of cluster #s ###
seq_cluster_options <- function(df, ev_list, start_num, stop_num, seed){
  ########
  # evs: list of vectors of eigenvector numbers (e.g., list(c(1), c(1,2), c(1,2,3)))
  # start_num: minimum cluster # to get results for
  # stop_num: maximum cluster # to get results for
  ########
  clust_fits <- lapply(seq(start_num, stop_num), function(x){
    do.call(cbind, lapply(ev_list, clust2, clust_num = x, df = df, seed = seed))})
  out_df <- cbind(df, do.call(cbind, clust_fits))
  return(out_df)
}
seq_hcluster_options <- function(df, dist_mat, start_num, stop_num, seed){
  clust_fits <- do.call(cbind, lapply(seq(start_num, stop_num), hclust_func, dist_mat = dist_mat, seed=seed))
  out_df <- cbind(df, clust_fits)
  return(out_df)
}



kmeans_mb_centers <- function(cols, df, k, seed){
  #print(df[cols])
  set.seed(seed)
  cols <- as.character(cols)
  clust_res <- kmeans(df[cols], k)
  col_name <- paste0('clust', k, '_', 'EV', paste0(cols, collapse = ''))
  colnames(clust_res$centers) <- paste0(col_name, '_center_', colnames(clust_res$centers))
  join_df <- as.data.frame(clust_res$centers) %>%
    tibble::rownames_to_column(col_name) %>%
    mutate_all(as.double)
  clust_res_df <- data.frame(clust_res$cluster)
  colnames(clust_res_df) <- col_name
  clust_res_df <- clust_res_df %>%
    left_join(join_df, by = col_name)
  return(clust_res_df)
}

seq_cluster_options_centers <- function(df, ev_list, start_num, stop_num, seed){
  #calc cluster fits for all combinations of evs and cluster #s to sequence through
  clust_fits <- lapply(seq(start_num, stop_num), function(x){
    do.call(cbind, lapply(ev_list, kmeans_mb_centers, k = x, df = df, seed = seed))})
  out_df <- cbind(df, do.call(cbind, clust_fits))
  return(out_df)
}







