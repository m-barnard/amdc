library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)

setwd('t2_AMDC_rev1/')
res_files <-paste0(seq(1, 20), '.rds')
all_res <- lapply( paste0('res_stability/',res_files), read_rds)
all_df <- do.call(rbind, all_res) 
all_df <- all_df %>%
  mutate(clust_num = str_sub(best_clust, 6, 6))
ggplot(all_df, aes(x = as.factor(clust_num))) +
  geom_bar() #mostly 7-9 clusters, most often 9 clusters

clust_nums <- all_df %>%
  dplyr::select(seed, clust_num) %>%
  unique()
sum(clust_nums$clust_num == 6)
491/500

comp_df <- read_csv('amdc_clust_res.csv', show_col_types = FALSE)
all_df <- all_df %>%
  left_join(comp_df, by = 'uniq_id')

agg_df1 <- all_df %>%
  group_by(uniq_id, clust8_EV1234) %>%
  summarise(jac = mean(jac)) %>%
  ungroup()
mean(agg_df1$jac) #0.7932195 is the overall score across all clusters, pretty good!
median(agg_df1$jac) #0.875684 is the median of this, again pretty good!!

#generally very good for the home and work cluster, pretty good for all other clusters except 2
agg_df2_med <- agg_df1 %>%
  group_by(clust8_EV1234) %>%
  summarise(jac = median(jac))
agg_df2_mean <- agg_df1 %>%
  group_by(clust8_EV1234) %>%
  summarise(jac = mean(jac))
ggplot(agg_df1, aes(x = jac, fill = as.factor(clust8_EV1234))) +
  geom_histogram() +
  facet_wrap(~clust8_EV1234, scales = 'free_y')









