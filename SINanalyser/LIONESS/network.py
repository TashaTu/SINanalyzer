#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 16:34:17 2024

@author: tashatu
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy.sparse
from scipy.sparse import csr_matrix 

def network_construction(data_df,mean,std,pvalue=0.05,sample_list=None,outdir=None,output_format="npz"):
    print('Step 3 : network construction')
    if not isinstance(data_df, pd.DataFrame):
        raise ValueError("`data_df` must be a Pandas DataFrame with genes as rows and samples as columns.")
    if not isinstance(mean, (int, float)):
        raise ValueError("`mean` must be a numeric value (int or float).")
    if not isinstance(std, (int, float)):
        raise ValueError("`standard` must be a numeric value (int or float).")
    if output_format not in ["npz", "edge_list", "edge_list_score", "edge_list_zscore"]:
        raise ValueError("`output_format` must be one of 'npz', 'edge_list', 'edge_list_score', or 'edge_list_zscore'.")
    if sample_list is not None and not all(sample in data_df.columns for sample in sample_list):
        raise ValueError("All samples in `sample_list` must exist in `data_df.columns`.")
        
    if not os.path.exists(outdir): os.makedirs(outdir) # Check if the directory exists
    zscore = stats.norm.isf(pvalue/2)
    # Data
    data_gene = list(data_df.index)
    data_sample = list(data_df.columns)
    data_val = data_df.values
    data_gene_l = len(data_gene)
    data_sample_l = len(data_sample)
    if sample_list == None: sample_list = data_sample
    print('-> Data preprocessing ... Done')
    
    if output_format == "npz":
        # Gene index mapping table
        with open(f'{outdir}/Gene_index_mapping_table.txt','w') as gene_file :
            for idx,val in enumerate(data_gene) :
                gene_file.write(f'{idx}\t{val}\n')
        print('-> Gene index mapping table ... Done')
    
    # Aggregate network construction
    agg = np.corrcoef(data_val)
    tri = np.full(agg.shape,False)
    tri[np.triu_indices(data_gene_l,1)] = True
    agg = agg[tri]
    agg[np.isnan(agg)] = 0
    print('-> Aggregate network construction ... Done')
    
    # Network contruction
    edge = np.full((data_gene_l,data_gene_l),False)
    for idx,sample in enumerate(data_sample):
        if sample in sample_list:
            edge_score = np.delete(data_val,idx,1)
            edge_score = np.corrcoef(edge_score)
            edge_score = edge_score[tri]
            edge_score[np.isnan(edge_score)] = 0
            edge_score =  data_sample_l * (agg - edge_score) + edge_score
            edge_score_z = np.abs((edge_score-mean)/std)
            edge.fill(False)
            edge[tri] = edge_score_z >= zscore
            rows, cols = np.where(edge)
    
            # Handle output formats
            if output_format == "npz":

                sparse_s = csr_matrix(edge)
                scipy.sparse.save_npz(f'{outdir}/{sample}.npz', sparse_s)
     
            elif output_format == "edge_list":
                with open(f'{outdir}/{sample}_edge_list.txt', 'w') as f:
                    for i, j in zip(rows, cols):
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\n')

            elif output_format == "edge_list_score":
                with open(f'{outdir}/{sample}_edge_list_score.txt', 'w') as f:
                    for i, j in zip(rows, cols):
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\t{edge_score[i, j]:.6f}\n')

            elif output_format == "edge_list_zscore":
                with open(f'{outdir}/{sample}_edge_list_zscore.txt', 'w') as f:
                    for i, j in zip(rows, cols):
                        f.write(f'{data_gene[i]}\t{data_gene[j]}\t{edge_score_z[i, j]:.6f}\n')
    print('-> Single-sample network construction ... Done')