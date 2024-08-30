################################
# Functions for creating the 'lasagna' plots
###############################

grab_state <- function(i, df){
  return(substr(df$seq, i, i))
}

create_lasagna_df <- function(col_name, df){
  seq_len <- nchar(df$seq[1])
  time_cols <- as.data.frame(sapply(seq(1, seq_len), grab_state, df = df))
  colnames(time_cols) <- as.character(seq(1, seq_len))
  seq_df <- df %>%
    select(uniq_id, sym(col_name)) %>%
    bind_cols(time_cols) %>%
    pivot_longer(as.character(seq(1, seq_len)), names_to ='Time', values_to = 'State') %>%
    mutate(Time = as.numeric(Time)) %>%
    mutate(O = ifelse(State == 'O',1,0), W = ifelse(State == 'W',1,0), H = ifelse(State == 'H',1,0),
           `T` = ifelse(State == 'T',1,0)) %>%
    select(-State) %>%
    pivot_longer(c('O', 'T', 'W', 'H'), names_to = 'State', values_to = 'Value')
  return(seq_df)
}



plot_lasagna <- function(seq_df, col_name,cluster){
  clust1 <- seq_df %>%
    filter(!!sym(col_name) == cluster) %>%
    mutate(int_var_w = ifelse((State == 'W') & (Value == 1),1,0)) %>%
    mutate(int_var_h = ifelse((State == 'H') & (Value == 1),1,0)) %>%
    mutate(int_var_t = ifelse((State == 'T') & (Value == 1),1,0)) %>%
    mutate(int_var_o = ifelse((State == 'O') & (Value == 1),1,0)) %>%
    group_by(uniq_id) %>%
    mutate(arrange_var = ifelse(sum(int_var_w) >80, 0, ifelse(sum(int_var_h) > 100, 1, ifelse(sum(int_var_t) > 200, 2, ifelse(sum(int_var_o) > 200,3,4))))) %>%
    arrange(arrange_var)
  
  join_df <- data.frame('uniq_id' = unique(clust1$uniq_id), 'num' = seq(0, length(unique(clust1$uniq_id))-1))
  clust1 <- clust1 %>%
    left_join(join_df) %>%
    mutate(Val_new = ifelse(Value == 0, NA,Value+num))

  ggplot(clust1) +
    geom_raster(aes(x = Time, y = Val_new, fill = State)) +
    theme_test() +
    scale_fill_manual(values = c("#80B1D3","#FDB462","#B3DE69","#FCCDE5")) +
    scale_x_continuous(breaks = c(0,72,144, 216,288), labels = c('12am', '6am', '12pm', '6pm', '12am')) +
    theme(axis.title = element_blank()) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), strip.background = element_rect(colour="white", fill="white"), panel.background = element_rect(colour="white", fill="white"), strip.text = element_text(face = "bold"), legend.position = 'bottom')
}



plot_lasagna_week <- function(seq_df, col_name,cluster){
  clust1 <- seq_df %>%
    filter(!!sym(col_name) == cluster) %>%
    mutate(int_var_w = ifelse((State == 'W') & (Value == 1),1,0)) %>%
    mutate(int_var_h = ifelse((State == 'H') & (Value == 1),1,0)) %>%
    mutate(int_var_t = ifelse((State == 'T') & (Value == 1),1,0)) %>%
    mutate(int_var_o = ifelse((State == 'O') & (Value == 1),1,0)) %>%
    group_by(uniq_id) %>%
    mutate(arrange_var = ifelse(sum(int_var_w) >80, 0, ifelse(sum(int_var_h) > 100, 1, ifelse(sum(int_var_t) > 200, 2, ifelse(sum(int_var_o) > 200,3,4))))) %>%
    arrange(arrange_var)
  
  join_df <- data.frame('uniq_id' = unique(clust1$uniq_id), 'num' = seq(0, length(unique(clust1$uniq_id))-1))
  clust1 <- clust1 %>%
    left_join(join_df) %>%
    mutate(Val_new = ifelse(Value == 0, NA,Value+num))
  
  ggplot(clust1) +
    geom_raster(aes(x = Time, y = Val_new, fill = State)) +
    theme_test() +
    geom_vline(xintercept = seq(289, 1440, 288), color = 'black')+
    scale_fill_manual(values = c("#80B1D3","#FDB462","#B3DE69","#FCCDE5")) +
    theme(axis.title = element_blank()) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), strip.background = element_rect(colour="white", fill="white"), panel.background = element_rect(colour="white", fill="white"), strip.text = element_text(face = "bold"), legend.position = 'bottom')
}


plot_lasagna_alpha <- function(seq_df, col_name,cluster){
  clust1 <- seq_df %>%
    filter(!!sym(col_name) == cluster) %>%
    mutate(int_var_w = ifelse((State == 'W') & (Value == 1),1,0)) %>%
    mutate(int_var_h = ifelse((State == 'H') & (Value == 1),1,0)) %>%
    mutate(int_var_t = ifelse((State == 'T') & (Value == 1),1,0)) %>%
    mutate(int_var_o = ifelse((State == 'O') & (Value == 1),1,0)) %>%
    group_by(uniq_id) %>%
    mutate(arrange_var = ifelse(sum(int_var_w) >80, 0, ifelse(sum(int_var_h) > 100, 1, ifelse(sum(int_var_t) > 200, 2, ifelse(sum(int_var_o) > 200,3,4))))) %>%
    arrange(diff_clust, arrange_var)
  
  join_df <- data.frame('uniq_id' = unique(clust1$uniq_id), 'num' = seq(0, length(unique(clust1$uniq_id))-1))
  clust1 <- clust1 %>%
    left_join(join_df) %>%
    mutate(Val_new = ifelse(Value == 0, NA,Value+num))
  
  ggplot(clust1) +
    geom_raster(aes(x = Time, y = Val_new, fill = State, alpha = diff_clust)) +
    theme_test() +
    scale_alpha_manual(values = c(1, 0.6)) +
    scale_fill_manual(values = c("#80B1D3","#FDB462","#B3DE69","#FCCDE5")) +
    scale_x_continuous(breaks = c(0,72,144, 216,288), labels = c('12am', '6am', '12pm', '6pm', '12am')) +
    theme(axis.title = element_blank()) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), strip.background = element_rect(colour="white", fill="white"), panel.background = element_rect(colour="white", fill="white"), strip.text = element_text(face = "bold"), legend.position = 'bottom')
}