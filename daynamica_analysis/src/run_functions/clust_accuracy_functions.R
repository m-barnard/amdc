#### Functions to assess simulated data method accuracy ###

calc_sum <- function(j, i, table){
  pull <- vector('numeric', length(i))
  for(h in i){
    pull[h] <- table[h, j[h]]
  }
  return(sum(pull))
}


calc_accuracy <- function(table, num_clusts){
  ####################
  # Calculate accuracy metric given a true and estimated set of clusters
  # table: table that compares the true and estimated clusters
  # num_clusts: the true # of clusters
  #####################
  index <- seq(1, num_clusts)
  max_vec <- apply(table,1, function(x){which(x == max(x))})
  
  #if the clustering method actually does a good job this if statement will be true
  if((length(unique(max_vec)) == length(max_vec)) & (typeof(max_vec) != 'list')){
    raw <- calc_sum(max_vec, index, table)/sum(table)
  }else{ #if clustering method bad, find the highest proportion/score you can give them though looking at possible permutations
    perm <- permn(seq(1,num_clusts))
    raw_best <- max(sapply(perm, calc_sum, i = index,table = table))
    raw <- raw_best/sum(table)
  }
  return(c('raw_prop' = raw))
}


get_mult_acc <- function(x_df){
  nt_tab <- table(as.factor(x_df$true_clust), as.factor(x_df$best_clust))
  prop_correct <- calc_accuracy(nt_tab, length(unique(x_df$true_clust)))
  return(c('best_num' = unique(x_df$best_num), 'seed' = unique(x_df$seed), prop_correct))
}


