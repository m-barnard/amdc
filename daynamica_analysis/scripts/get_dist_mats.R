library(readr)

df <- read_csv('clean_data/day_seqs.csv')
out <- adist(df$seqs)
write_rds(out, 'clean_data/day_distmat.rds')

week_df <- read_csv('clean_data/clean_week_seqs.csv')
out <- adist(week_df$seqs)
write_rds(out, 'clean_data/clean_week_distmat.rds')

