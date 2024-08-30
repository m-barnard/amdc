import numpy as np
import pandas as pd
import os
import itertools
from itertools import repeat
import sys
sys.path.insert(1,  f'{os.path.abspath(os.path.join(os.getcwd() ,"../.."))}/src/run_functions')
import run_nTreeClus_one_df
from run_nTreeClus_one_df import get_clust_res
import multiprocessing


######################
# Run nTreeClus for all datasets in a simulation scenario
# NOTE: simulation inputs are mainly created through reading in folder names for the generated simulated sequences
# so the inputs may require some tweaking if all sequences are not simulated before this
# i: indexes the input_grid matrix to identify the simulation scenario to run
# output: all data is saved as csv within results/order_{}/ and named by simulation scenario
######################

#keep this setting if running code on a slurm cluster, otherwise set i manually 1-36
i = int(sys.argv[1])
path = os.path.abspath(os.path.join(os.getcwd() ,"../.."))

logging.info(f'start {i}')
outer_folders = sorted(os.listdir(f'{path}/simulated_sequences/'))
folders = sorted(os.listdir(f'{path}/simulated_sequences/order_1/'))
input_grid = list(itertools.product(outer_folders,folders))
input_grid = pd.DataFrame(input_grid,columns=("o_fold","folder_name"))
input_grid = pd.concat([input_grid, pd.DataFrame([6]*len(input_grid.index), columns = ['max']),pd.DataFrame([3]*len(input_grid.index), columns = ['true'])], axis = 1)
input_grid.loc[input_grid['folder_name'] == 'state_high', 'max'] = 7
input_grid.loc[input_grid['folder_name'] == 'state_low', 'max'] = 5
input_grid.loc[input_grid['folder_name'] == 'state_high', 'true'] = 4
input_grid.loc[input_grid['folder_name'] == 'state_low', 'true'] = 2

inputs = input_grid.iloc[i]
   
#run ntreeClus for all datasets in a given simulation scenario
pool = multiprocessing.Pool(24)
res = pool.starmap(get_clust_res, zip(range(0, 500), repeat(path), repeat(inputs)))
final_df = pd.concat(res, axis = 0)
final_df.to_csv(f'{path}/results/{inputs["o_fold"]}/{inputs["folder_name"]}_ntrees_clust.csv', index=False)
    
