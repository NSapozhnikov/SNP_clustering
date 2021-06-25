import numpy as np
import pandas as pd
import sys
import argparse  

def __main__():
    return
__main__()

snp_file_path = 'chr1_snplist.snplist'
# kmeans_np = np.loadtxt('chr1_kmeans.clustered', dtype=str) 
# dbscan_np = np.loadtxt('chr1_dbscan.clustered', dtype=str) 
# hdbscan_np = np.loadtxt('chr1_hdbscan.clustered', dtype=str)
# spectral_np = np.loadtxt('chr1_spectral.clustered', dtype=str)

kmeans_df = pd.read_table('chr1_kmeans.clustered', sep=' ', header=None)
dbscan_df = pd.read_table('chr1_dbscan.clustered', sep=' ', header=None)
hdbscan_df = pd.read_table('chr1_hdbscan.clustered', sep=' ', header=None)
spectral_df = pd.read_table('chr1_spectral.clustered', sep=' ', header=None)

num_clust_kmeans = int(kmeans_df.max(numeric_only=True) + 1)
num_clust_dbscan = int(dbscan_df.max(numeric_only=True) + 2)
num_clust_hdbscan = int(hdbscan_df.max(numeric_only=True) + 2)
num_clust_spectral = int(spectral_df.max(numeric_only=True) + 1)

# method = input('File 1 to compare:')

# num_clust_1 = dbscan_df.max(numeric_only=True) + (1 if method == ('kmeans' or 'specral') else 2)

with open (snp_file_path, 'r') as snp_file: 
    snp_list = np.loadtxt(snp_file, dtype=str)
    

with open('kmeans_spectral.comparison', 'w') as out_file:
    for i in range(num_clust_kmeans):
        clu = kmeans_df.loc[kmeans_df[1] == i]#.reset_index(drop=True)
        for j in range(num_clust_spectral):        
            clu1 = spectral_df.loc[spectral_df[1] == j]#.reset_index(drop=True)
            clu_merged = clu.merge(clu1, left_on=0, right_on=0, how='left')
            if not clu_merged['1_y'].isnull().values.any():
                out_file.write(clu_merged.to_string())
                out_file.write('\n')
                       

