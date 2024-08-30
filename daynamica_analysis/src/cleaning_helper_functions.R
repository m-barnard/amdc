split_mult_days <- function(df){
  ############################
  # Split calendar items that go through multiple days
  ############################
  id_df <- df %>%
    mutate(mult_day = ifelse(date(start_posix) != date(end_posix),1,0))
  single_df <- id_df %>%
    filter(mult_day == 0) %>%
    dplyr::select(-mult_day)
  mult_df <- id_df %>%
    filter(mult_day == 1) %>%
    dplyr::select(-mult_day) %>%
    mutate(num_days = as.double(date(end_posix) -date(start_posix) + 1))
  mult_df_exp <- mult_df %>%
    mutate(rn = row_number()) %>%
    rowwise() %>%
    slice(rep(1, num_days)) %>%
    group_by(rn) %>%
    #the 60*60*24 adds a day, floor gets 00:00:00
    mutate(start_posix = c(start_posix[1], seq(floor_date((start_posix[1] + 60*60*24), unit = "days"), 
                                               floor_date(end_posix[1], unit = "days"),
                                               by = "1 day"))) %>%
    #subtract one to get the 23:59:59
    mutate(end_posix = c(seq(floor_date(start_posix[1] + 60*60*24, unit = "days"), 
                             floor_date(end_posix[1], unit = "days"),
                             by = "1 day") -1, end_posix[1])) %>%
    ungroup() %>%
    dplyr::select(-rn, -num_days)
  final_df <- bind_rows(single_df, mult_df_exp) %>% arrange(user_id, start_posix, end_posix, cal_item_id)
  return(final_df)
}


time_to_sequence <- function(row, segments, segment_midpoints){
  ######################
  # Convert the calendar items to sequence
  ####################
  day_sequence <- rep(NA, segments)
  day_sequence[which(segment_midpoints >= as.numeric(row['start_minute_total']) & segment_midpoints <= as.numeric(row['end_minute_total']))] <- row['subtype_trunc']
  new_seq <- day_sequence
  return(new_seq)
}

day_to_sequence <- function(day_data, segments, segment_midpoints){
  ##############
  # Convert a day of calendar items into a sequence
  ##############
  test2 <-apply(day_data, 1, time_to_sequence, segments = segments, segment_midpoints = segment_midpoints)
  vec_list <- lapply(seq_len(ncol(test2)), function(i) test2[,i])  
  final_vec <- coalesce(!!!vec_list)
  final_vec[is.na(final_vec)] <- "X"
  seq <- paste0(final_vec, collapse = '')
  return(seq)
}

data_to_sequence <- function(data, bandwidth, conversion_codes){
  #####################
  # Get set of sequences for a given dataset
  # bandwidth: number of minutes each letter in the sequence should represent
  #####################
  segments = round(1440/bandwidth)
  segment_midpoints <- 1:segments*bandwidth - (bandwidth/2)
  data <- data %>% 
    mutate(start_minute_total = hour(start_posix)*60 + minute(start_posix),
           end_minute_total = hour(end_posix)*60 + minute(end_posix),
           subtype = as.character(subtype),
           subtype_trunc = conversion_codes[subtype],
           uniq_id = paste(user_id, '_',as.character(date(start_posix))))
  
  unique_id_dates <- data %>%
    mutate(start_date = date(start_posix)) %>%
    select(user_id, start_date, uniq_id) %>%
    unique()
  
  list_data <- split(data, f = data$uniq_id)
  final_seq <- lapply(list_data, day_to_sequence, segments = segments, segment_midpoints = segment_midpoints)
  final_df <- data.frame(Reduce(rbind, final_seq)) %>%
    rename(seq = `Reduce.rbind..final_seq.`) %>%
    mutate(uniq_id = names(final_seq)) %>%
    full_join(unique_id_dates, by = 'uniq_id')
  
  return(final_df)
}