import numpy as np
import datetime
from sklearn.cluster import  DBSCAN, KMeans, SpectralClustering 
from sklearn.metrics.cluster import *
import hdbscan
import pandas as pd
import sys
import argparse    
from getpass import getpass
import re
import configparser
from joblib import Memory

def dbscan_clustering(diss_matrix, dbscan_shil_score_list, dbscan_CH_list):
    print('\nPerforming clustering with DBSCAN...')
    db = DBSCAN(eps=0.75, min_samples=6, metric='precomputed', n_jobs=10).fit(diss_matrix)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    ### Scoring
    try:
        dbscan_shil_score = silhouette_score(diss_matrix, labels, metric='precomputed') 
        print('Silhouette score: ', dbscan_shil_score)
        dbscan_shil_score_list.append(dbscan_shil_score)
    except:
        print('Failed to compute Silhouette score.')
        dbscan_shil_score_list.append('-')
    try:
        dbscan_CH = calinski_harabasz_score(diss_matrix, labels)
        print('Calinski-Harabasz score: ', dbscan_CH)
        dbscan_CH_list.append(dbscan_CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
        dbscan_CH_list.append('-')
    return dbscan_shil_score_list, dbscan_CH_list
        
def hdbscan_clustering(diss_matrix, hdbscan_shil_score_list, hdbscan_CH_list):
    print('\nPerforming clustering with HDBSCAN...')
    mem = Memory(cachedir='./cachedir')
    clusterer = hdbscan.HDBSCAN(min_cluster_size=8, min_samples=6, 
                                metric='precomputed', cluster_selection_method='leaf',
                                core_dist_n_jobs=-1, memory=mem).fit(diss_matrix)
    labels = clusterer.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    ### Scoring
    try:
        hdbscan_shil_score = silhouette_score(diss_matrix, labels, metric='precomputed') 
        print('Silhouette score: ', hdbscan_shil_score)
        hdbscan_shil_score_list.append(hdbscan_shil_score)
    except:
        print('Failed to compute Silhouette score.')
        hdbscan_shil_score_list.append('-')
    try:
        hdbscan_CH = calinski_harabasz_score(diss_matrix, labels)
        print('Calinski-Harabasz score: ', hdbscan_CH)
        hdbscan_CH_list.append(hdbscan_CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
        hdbscan_CH_list.append('-')
    mem.clear(warn=False)
    return hdbscan_shil_score_list, hdbscan_CH_list

def kmeans_clustering(corr_matrix, k_shil_score_list, k_CH_list):
    print('\nPerforming clustering with K-means...')
    kmeans = KMeans(n_clusters=100, n_init=10, max_iter=100).fit(corr_matrix)
    labels = kmeans.predict(corr_matrix)
    n_clusters_ = len(set(labels))
    print('Chosen number of clusters: %d' % n_clusters_)
    ### Scoring
    try:
        k_shil_score = silhouette_score(corr_matrix, labels) 
        print('Silhouette score: ', k_shil_score)
        k_shil_score_list.append(k_shil_score)
    except:
        print('Failed to compute Silhouette score.')
        k_shil_score_list.append('-')
    try:
        k_CH = calinski_harabasz_score(corr_matrix, labels)
        print('Calinski-Harabasz score: ', k_CH)
        k_CH_list.append(k_CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
        k_CH_list.append('-')
    return k_shil_score_list, k_CH_list
    
def spectral_clustering(corr_matrix, spectral_shil_score_list, spectral_CH_list):
    print('\nPerforming clustering with SpectralClustering...')
    clustering = SpectralClustering(n_clusters=100, affinity='precomputed', n_jobs=-1).fit(corr_matrix)
    labels = clustering.labels_
    n_clusters_ = len(set(labels))
    print('Chosen number of clusters: %d' % n_clusters_)
    ### Scoring
    try:
        spectral_shil_score = silhouette_score(corr_matrix, labels)
        print('Silhouette score: ', spectral_shil_score)
        spectral_shil_score_list.append(spectral_shil_score)
    except:
        print('Failed to compute Silhouette score.')
        spectral_shil_score_list.append('-')
    try:
        spectral_CH = calinski_harabasz_score(corr_matrix, labels)
        print('Calinski-Harabasz score: ', spectral_CH)
        spectral_CH_list.append(spectral_CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
        spectral_CH_list.append('-')
    return spectral_shil_score_list, spectral_CH_list
  
def main_():
    k_shil_score_list, k_CH_list, dbscan_shil_score_list = [], [], []
    dbscan_CH_list, hdbscan_shil_score_list, hdbscan_CH_list = [], [], []  
    spectral_shil_score_list, spectral_CH_list, chr_list = [], [], []
    for num_chr in range(1,23):
        snp_file_path = f"/mnt/wd/nsap/old/input_files/chr{num_chr}_old.snplist"
        matrix_file_path = f"/mnt/wd/nsap/old/input_files/chr{num_chr}_squared_old.ld"
        try:
            with open(snp_file_path, 'r') as snp_file:
              print('\nExtracting SNP list data from file...')
              snp_list = np.loadtxt(snp_file, dtype=str)        
              print(f"The SNP list for {str(num_chr)} chromosome is: \n{snp_list}")
              print(f"With length of: {len(snp_list)}")
        except FileNotFoundError:
            sys.exit('No such file.')      
        try:
            with open(matrix_file_path, 'r') as corr_file:
                print('Extracting matrix data from file...')
                corr_matrix = np.fromfile(corr_file, sep=' ')
        except FileNotFoundError:
            sys.exit('No such file.')
        np.nan_to_num(corr_matrix, copy=False)  
        print('Reshaping an array...')
        diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
        diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
        np.fill_diagonal(diss_matrix, 0)
        corr_matrix = corr_matrix.reshape(len(snp_list), len(snp_list))
        np.fill_diagonal(corr_matrix, 1)
        ### Clustering
        ## k-means
        kmeans_clustering(corr_matrix, k_shil_score_list, k_CH_list)
        ## dbscan
        dbscan_clustering(diss_matrix, dbscan_shil_score_list, dbscan_CH_list)
        ## hdbscan
        hdbscan_clustering(diss_matrix, hdbscan_shil_score_list, hdbscan_CH_list)
        ## spectral
        spectral_clustering(corr_matrix, spectral_shil_score_list, spectral_CH_list)
        ### Add index
        chr_list.append(num_chr)
    ### Merging results
    pd_df = pd.DataFrame(data={'chr': chr_list,
                           'dbscan_shil_score': dbscan_shil_score_list,
                           'dbscan_CH': dbscan_CH_list,
                           'hdbscan_shil_score': hdbscan_shil_score_list,
                           'hdbscan_CH': hdbscan_CH_list,
                           'k_shil_score': k_shil_score_list,
                           'k_CH': k_CH_list,
                           'spectral_shil_score': spectral_shil_score_list,
                           'spectral_CH': spectral_CH_list})
    pd_df.to_excel('scores_old.xlsx', index=False)
main_()

        
        
        
        
