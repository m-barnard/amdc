grab_state <- function(i, df){
  return(substr(df$seq, i, i))
}


create_lasagna_df <- function(df, states){
  col_name <- 'cluster'
  seq_len <- nchar(df$seq[1])
  
  count <- df %>%
    group_by(!!sym(col_name)) %>%
    summarize(n = n())

  set.seed(5)
  df <- df %>%
    group_by(!!sym(col_name)) %>%
    sample_n(size = min(count$n), replace = FALSE)
  
  time_cols <- as.data.frame(sapply(seq(1, seq_len), grab_state, df = df))
  colnames(time_cols) <- as.character(seq(1, seq_len))
  
  seq_df <- df %>%
    select(uniq_id, sym(col_name), sim) %>%
    bind_cols(time_cols) %>%
    pivot_longer(as.character(seq(1, seq_len)), names_to ='Time', values_to = 'State') %>%
    mutate(Time = as.numeric(Time))
  for(k in states){
    seq_df <- cbind(seq_df, ifelse(seq_df$State == k,1,0))
    colnames(seq_df) <- c(colnames(seq_df)[1:length(colnames(seq_df)) -1], k)
  }
  
  seq_df <- seq_df %>%
    dplyr::select(-State) %>%
    pivot_longer(states, names_to = 'State', values_to = 'Value') %>%
    arrange(cluster,uniq_id, Time)
  return(seq_df)
}



plot_lasagna <- function(seq_df, states){
  col_name <- 'cluster'
  df_list <- list()
  seq_df <- seq_df %>%
    # mutate(sim_type =  paste(sim, type, '\n')) %>%
    mutate(sim_type = paste0(sim, ' Overlap')) %>%
    mutate(sim_type = factor(sim_type, levels = c('Low Overlap', 'Medium Overlap', 'High Overlap')))
  #mutate(sim_type = factor(sim_type, levels = paste(overlap_names, type_names, '\n')))
  index <- seq_df %>%
    dplyr::select(!!sym(col_name), sim_type) %>%
    unique()
  for(i in seq(1, nrow(index))){
    k <- index[i, ]$cluster
    j <- index[i, ]$sim_type
    clust1 <- seq_df %>%
      filter(!!sym(col_name) == k & sim_type == j)
    join_df <- data.frame('uniq_id' = unique(clust1$uniq_id), 'num' = seq(0, length(unique(clust1$uniq_id))-1))
    clust1 <- clust1 %>%
      left_join(join_df) %>% 
      arrange(uniq_id) %>%
      mutate(Val_new = ifelse(Value == 0, NA,Value+num)) 
    df_list[[i]] <- clust1
  }
  
  clust1 <- do.call(rbind, df_list)
  if(length(states) ==4){
    colors <- c("#80B1D3","#FDB462","#B3DE69","#FCCDE5")
  } else if(length(states) ==3){
    colors <- c("#80B1D3","#FDB462","#B3DE69")
  }
  
  ggplot(clust1) +
    geom_raster(aes(x = Time, y = Val_new, fill = State), position = 'identity') +
    theme_test() +
    scale_fill_manual(values = colors) +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line.y = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), strip.background = element_rect(colour="white", fill="white"), panel.background = element_rect(colour="white", fill="white"), strip.text = element_text(face = "bold"), legend.position = 'bottom') +
    facet_grid(vars(sim_type),vars(cluster), scales = 'free', labeller = label_wrap_gen(width=10))
}
