import numpy as np
import pandas as pd
import os
import itertools
from itertools import repeat
import sys
from ntrees_edit import nTreeClus

def get_clust_res(df_num, path, inputs):
    ''' Inputs: 
            df_num: simulated dataset #
            path: path to simulations/ folder 
            inputs: a dictionary of the following form
            {'o_fold':{outer folder name}, 'folder_name':{inner folder name}, 'max': {max # of clusters}, 'true': {true # of clusters}}
    '''
    #find best # of clusters
    df = pd.read_csv(f'{path}/simulated_sequences/{inputs["o_fold"]}/{inputs["folder_name"]}/{df_num}.csv')
    sequences = df.loc[:,'seqs'].tolist()
    model     = nTreeClus(sequences, n = None, ntree=10, method = "All", max_clust = inputs["max"], verbose = 0)
    model.nTreeClus()
    out = model.output()
    ord_names = ['DT', 'RF', 'DT_p', 'RF_p']
    sil_list = [out['sil_DT'], out['sil_RF'], out['sil_DT_p'], out['sil_RF_p']]
    best = ord_names[sil_list.index(max(sil_list))]
    best_num = out[f'C_{best}']
    
    #get best set clusters for the true# of clusters
    model2     = nTreeClus(sequences, n = None, ntree=10, method = "All", C = inputs["true"],verbose = 0)
    model2.nTreeClus()
    out2 = model2.output()
    ord_names2 = ['DT', 'RF', 'DT_p', 'RF_p']
    sil_list2 = [out2['sil_DT'], out2['sil_RF'], out2['sil_DT_p'], out2['sil_RF_p']]
    best2 = ord_names2[sil_list2.index(max(sil_list2))] # if more than one max, always goes with the first (based on ordering this will be simpler submethod)
    best_clust = out2[f'labels_{best2}']

    best_df = pd.DataFrame({'best_num': [best_num]*len(best_clust), 'best_clust': best_clust, 'true_clust': df.loc[:,'cluster'].tolist(), 'seed': [df_num]*len(best_clust)})
    return(best_df)
    return(print(os.getcwd()))