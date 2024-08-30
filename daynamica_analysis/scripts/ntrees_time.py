import numpy as np
import pandas as pd
import os
from itertools import repeat
import sys
sys.path.insert(1, f'{os.path.abspath(os.path.join(os.getcwd() ,"../"))}/src/run_functions')
import ntrees_edit
from ntrees_edit import nTreeClus
import itertools

################
# Get computational time it takes to run nTreeClus
# Run each method 25 times and take mean of the 25 times in agg_time_results.R
################

path = os.path.abspath(os.path.join(os.getcwd() ,"../"))
prep_time = []
all_time = []
for k in range(0,25):
    df = pd.read_csv(f'{path}/clean_data/clean_day_seqs.csv')
    sequences = df.loc[:,'seqs'].tolist()
    model     = nTreeClus(sequences, n = None, ntree=10, method = "All", max_clust = 9, verbose = 0)
    model.nTreeClus()
    out = model.output()
    prep_k = out['running_timeSegmentation']
    all_k = prep_k + out['running_timeDT'] + out['running_timeDT_p'] + out['running_timeRF'] + out['running_timeRF_p']
    prep_time += [prep_k]
    all_time += [all_k]
    
with open(f'{path}/results_time/ntrees_prep_day.txt', 'w') as f:
    for t in prep_time:
        f.write(f"{t}\n")
with open(f'{path}/results_time/ntrees_day.txt', 'w') as f:
    for t in all_time:
        f.write(f"{t}\n")

prep_time2 = []
all_time2 = []
for k in range(0, 25):
    df2 = pd.read_csv(f'{path}/clean_data/clean_week_seqs.csv')
    sequences2 = df2.loc[:,'seqs'].tolist()
    model2     = nTreeClus(sequences2, n = None, ntree=10, method = "All", max_clust = 9, verbose = 0)
    model2.nTreeClus()
    out2 = model2.output()
    prep_k = out2['running_timeSegmentation']
    all_k = prep_k + out2['running_timeDT'] + out2['running_timeDT_p'] + out2['running_timeRF'] + out2['running_timeRF_p']
    prep_time2 += [prep_k]
    all_time2 += [all_k]
    
with open(f'{path}/results_time/ntrees_prep_week.txt', 'w') as f:
    for t in prep_time2:
        f.write(f"{t}\n")
with open(f'{path}/results_time/ntrees_week.txt', 'w') as f:
    for t in all_time2:
        f.write(f"{t}\n")


    
