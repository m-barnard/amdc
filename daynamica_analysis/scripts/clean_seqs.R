library(dplyr)
library(readr)
library(tibble)
library(stringr)
library(lubridate)
source('src/cleaning_helper_functions.R')

### Clean day sequences ###
#clean calendar data
p_df <-read_csv("raw_data/processed_ucal_items_without_surveys_small.csv")
labels <- c("T", "T", "T",
            "W", "H", "O",
            "O", "T", "O",
            "O", "O", "T",
            "T", "T", "T",
            "W","T","X")
names(labels) <- c("OTHER","WALK", "CAR",
                   "WORK", "HOME", "SHOP",
                   "UNKNOWN_ACTIVITY", "IN_VEHICLE", "LEISURE_RECREATION", 
                   "EAT_OUT", "PERSONAL_BUSINESS","WAIT",
                   "BIKE", "BUS","UNKNOWN_TRAVEL_MODE",
                   "EDUCATION", "RAIL", "MISSING")
out <- data_to_sequence(p_df, 5, labels)


#remove missing data, and for individuals with more than 20 sequences, randomly select only 20
out <- out %>%
  rowwise() %>%
  mutate(num_states = length(unique(strsplit(seq, '')[[1]]))) %>%
  mutate(num_missing = sum(str_detect('X', strsplit(seq, '')[[1]]))) %>%
  ungroup() %>%
  mutate(hour_missing = num_missing*5/60) %>%
  mutate(missing_ind = as.factor(ifelse(hour_missing > 0, 1, 0))) %>%
  mutate(weekday = weekdays(start_date)) %>%
  mutate(day_type = ifelse(weekday == 'Saturday' | weekday == 'Sunday', 'weekend', 'M-F'))

df_prep <- out %>%
  filter(missing_ind ==0) %>%
  group_by(user_id) %>%
  mutate(num_seq_ind = n()) %>%
  ungroup()
df1 <- df_prep %>%
  filter(num_seq_ind <= 20)
#randomly sample 20 rows from each individual with more than 20 rows
set.seed(5)
df2 <- df_prep %>%
  filter(num_seq_ind > 20) %>%
  group_by(user_id) %>%
  sample_n(20)
df <- bind_rows(df1, df2) %>%
  rename(seqs = seq)

get_state_count <- function(x, state){
  vec <- unlist(strsplit(x, split = ""))
  return(sum(vec == state))
}

df$num_O <- sapply(df$seqs, get_state_count, state = 'O')
df$num_T <- sapply(df$seqs, get_state_count, state = 'T')
df$num_W <- sapply(df$seqs, get_state_count, state = 'W')

full_df <- df %>%
  #90% of sequence cannot be one of O, T, W
  filter(num_O < 260 & num_T < 260 & num_W < 260)
write_csv(full_df, 'clean_data/clean_day_seqs.csv')


### Clean week sequences ###
full_df_exp <- full_df %>%
  arrange(user_id, start_date) %>%
  group_by(user_id) %>%
  mutate(date_diff = c(NA, diff(start_date))) %>%
  mutate(date_lag = lead(date_diff)) %>%
  mutate(final_lag = coalesce(date_lag,date_diff))

fill_correct <- function(vec){
  final_vec <- vector(length = length(vec))
  fill = 1
  for(i in seq(1, length(vec))){
    if(vec[i] == 1){
      final_vec[i] <- fill
    }else{
      final_vec[i] <- fill
      fill <- fill + 1
    }
  }
  return(final_vec)
}

#count the number of consecutive days per individual
consec_days <- full_df_exp %>%
  filter(date_diff == 1 | date_lag == 1) %>%
  group_by(user_id) %>%
  mutate(count_consec = fill_correct(final_lag)) %>%
  group_by(user_id, as.factor(count_consec)) %>%
  mutate(consec_length = n()) %>%
  filter(consec_length > 4) 

fill_dates <- function(vec){
  final_vec <- vector(length = length(vec))
  fill = FALSE
  for(i in seq(1, length(vec))){
    if(vec[i] != 'Monday'){
      final_vec[i] <- fill
    }else{
      fill = TRUE
      final_vec[i] <- fill
    }
  }
  return(final_vec)
}


#each of the following three datasets gets a week sequence for a specific type of data
#week_sub where the consecutive Monday-Friday is broken with a weekend
#week_full_standard where there is one consecutive Monday-Friday
#week_more where there is two consective Monday-Friday sequences per individual
week_sub <- consec_days %>%
  filter(consec_length < 7) %>%
  group_by(user_id, as.factor(count_consec)) %>% #need to have count consec bc some people have two separate consec sequences > 4
  mutate(tot_week = sum(day_type == 'M-F')) %>%
  filter(tot_week == 5 & day_type != 'weekend') %>%
  mutate(weekday = factor(weekday, levels= c("Monday", 
                                             "Tuesday", "Wednesday", "Thursday", "Friday"))) %>%
  arrange(user_id, as.factor(count_consec), weekday) %>%
  summarize(seqs = paste0(seqs, collapse = "")) %>%
  ungroup() %>%
  dplyr::select(seqs, user_id)

week_full_standard <- consec_days %>%
  filter(consec_length > 6 & consec_length < 14) %>%
  filter(day_type != 'weekend') %>% 
  arrange(user_id, as.factor(count_consec), start_date) %>%
  group_by(user_id, as.factor(count_consec)) %>%
  slice(1:5) %>%
  mutate(weekday = factor(weekday, levels= c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))) %>%
  arrange(user_id, as.factor(count_consec), weekday) %>%
  summarize(seqs = paste0(seqs, collapse = "")) %>%
  ungroup() %>%
  dplyr::select(seqs, user_id)

week_more <- consec_days %>%
  filter(consec_length > 13) %>%
  group_by(user_id, as.factor(count_consec)) %>%
  mutate(vals = seq(1, unique(consec_length))) %>%
  mutate(sec = ifelse(vals < 8, 1, ifelse(vals < 15, 2, NA))) %>%
  filter(is.na(sec) == FALSE) %>%
  filter(day_type != 'weekend') %>% 
  group_by(user_id, as.factor(count_consec), sec) %>%
  slice(1:5) %>%
  mutate(weekday = factor(weekday, levels= c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))) %>%
  arrange(user_id, as.factor(count_consec), sec, weekday) %>%
  summarize(seqs = paste0(seqs, collapse = "")) %>%
  ungroup() %>%
  dplyr::select(seqs, user_id)



weeks_final <- bind_rows(week_sub, week_full_standard, week_more)
write_csv(weeks_final,'clean_data/clean_week_seqs.csv')

