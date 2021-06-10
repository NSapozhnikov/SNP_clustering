"""snp_clustering v0.5
Please check README.txt"""

import numpy as np
import datetime
from sklearn.cluster import AffinityPropagation, DBSCAN, KMeans, SpectralClustering
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import matplotlib.pyplot as plt
from itertools import cycle
from sklearn.metrics import pairwise_distances
import hdbscan
import seaborn as sns
from sklearn_extra.cluster import KMedoids
import pandas as pd
import sys
import argparse    
# """to do: implement data plotting"""

def hierarchy_clustering(diss_matrix_triangled, matrix_file_path):
    print('\nPerforming clustering with scipy.cluster.hierarchy...')
    print(f"Triangled dissimilarity matrix is: {diss_matrix_triangled}")
    print(f"With length of: {len(diss_matrix_triangled)}")
    hierarchy = linkage(diss_matrix_triangled, method='single')
    print('\nOutput matrix is: ')
    print(hierarchy)

def dbscan_clustering(diss_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with DBSCAN...')
    db = DBSCAN(eps=0.7, min_samples=6, metric='precomputed').fit(diss_matrix)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    out_file_path = matrix_file_path.rstrip('.ld')
    with open(f"{out_file_path}_dbscan.clustered", 'w') as fp:
        print('Writting to outfile...')
        for i in range(-1, n_clusters_):
            cluster_indices = np.where(labels == i)[0]
            for j in cluster_indices:
                fp.write(str(snp_list[j]))
                fp.write(' ')
                fp.write(str(i))
                fp.write('\n') 

def af_prop_clustering(corr_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with AffinityPropagation...')
    clustering = AffinityPropagation(damping=0.99, max_iter=200,
                                     convergence_iter=15, 
                                     affinity='precomputed', random_state=0).fit(corr_matrix)
    cluster_centers_indices = clustering.cluster_centers_indices_
    labels = clustering.labels_
    n_clusters_ = len(cluster_centers_indices)
    print('Estimated number of clusters: %d' % n_clusters_)
    out_file_path = matrix_file_path.rstrip('.ld')
    with open(f"{out_file_path}_afprop.clustered", 'w') as fp:
        print('Writting to outfile...')
        for i in range(n_clusters_):
            cluster_indeces = np.where(labels == i)[0]
            for j in cluster_indeces:
                fp.write(str(snp_list[j]))
                fp.write(' ')
                fp.write(str(i))
                fp.write('\n')

def hdbscan_clustering(corr_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with HDBSCAN...')
    clusterer = hdbscan.HDBSCAN(min_cluster_size=8, min_samples=2, metric='precomputed', cluster_selection_method='leaf').fit(corr_matrix)
    labels = clusterer.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    out_file_path = matrix_file_path.rstrip('.ld')
    with open(f"{out_file_path}_hdbscan.clustered", 'w') as fp:
        print('Writting to outfile...')
        for i in range(-1, n_clusters_):
            cluster_indeces = np.where(labels == i)[0]
            for j in cluster_indeces:
                fp.write(str(snp_list[j]))
                fp.write(' ')
                fp.write(str(i))
                fp.write('\n')

def kmeans_clustering(corr_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with K-means...')
    kmeans = KMeans(n_clusters=8, n_init=100).fit(corr_matrix)
    labels = kmeans.predict(corr_matrix)
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    out_file_path = matrix_file_path.rstrip('.ld')
    with open(f"{out_file_path}_kmeans.clustered", 'w') as fp:
        print('Writting to outfile...')
        for i in range(n_clusters_):
            cluster_indeces = np.where(labels == i)[0]
            for j in cluster_indeces:
                fp.write(str(snp_list[j]))
                fp.write(' ')
                fp.write(str(i))
                fp.write('\n')
    
def kmedoids_clustering(corr_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with K-medoids...')
    kmedoids = KMedoids(n_clusters=8, metric='correlation', method='pam', max_iter=1000).fit(corr_matrix)
    labels = kmedoids.predict(corr_matrix)
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    out_file_path = matrix_file_path.rstrip('.ld')
    with open(f"{out_file_path}_kmedoids.clustered", 'w') as fp:
        print('Writting to outfile...')
        for i in range(n_clusters_):
            cluster_indeces = np.where(labels == i)[0]
            for j in cluster_indeces:
                fp.write(str(snp_list[j]))
                fp.write(' ')
                fp.write(str(i))
                fp.write('\n')
                
def spectral_clustering(corr_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with SpectralClustering...')
    clustering = SpectralClustering(n_clusters=8, affinity='precomputed').fit(corr_matrix)
    labels = clustering.labels_
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    out_file_path = matrix_file_path.rstrip('.ld')
    with open(f"{out_file_path}_spectral.clustered", 'w') as fp:
        print('Writting to outfile...')
        for i in range(n_clusters_):
            cluster_indeces = np.where(labels == i)[0]
            for j in cluster_indeces:
                fp.write(str(snp_list[j]))
                fp.write(' ')
                fp.write(str(i))
                fp.write('\n') 

def awaiting_input():
    try:
        matrix_file_path = str(sys.argv[1])
        f = open(matrix_file_path, 'r')
        f.close()        
    except (IndexError, FileNotFoundError, ValueError ):
        sys.exit('Invalid matrix file path.')
    try:
        snp_file_path = str(sys.argv[2])
        f = open(snp_file_path, 'r')
        f.close()
    except (IndexError, FileNotFoundError, ValueError ):
        sys.exit('Invalid SNP list file path.')    
    try:
        methods = ['hierarchy', 'dbscan', 'affinity_propagation',
               'hdbscan', 'k_means', 'k_medoids', 'spectral']
        result = 'NOT OK'
        i = 0
        typed_input = sys.argv[3].lower().replace('-', '_')
        for met in methods:
            if typed_input != met:
                i+=1                
            elif typed_input == '':
                i+=7         
            else:
                print(f"\n{met} algorithm will be used...")
                time_start = datetime.datetime.now()
                print (f"\nTime of start - {time_start.isoformat(sep=' ', timespec='seconds')}")    
                try:
                    with open (snp_file_path, 'r') as snp_file: 
                        print('\nExtracting SNP list data from file...')
                        snp_list = np.loadtxt(snp_file, dtype=str)
                except FileNotFoundError:
                    sys.exit('No such file.')
                try:
                    with open(matrix_file_path, 'r') as corr_file:
                        print('Extracting matrix data from file...')
                        corr = np.fromfile(corr_file, sep=' ')
                except FileNotFoundError:
                    sys.exit('No such file.')
                print('Reshaping an array...')
                corr_matrix = corr.reshape(len(snp_list), len(snp_list))
                dissimilarity = 1 - abs(corr)
                diss_matrix = dissimilarity.reshape(len(snp_list), len(snp_list))
                diss_matrix_triangled = squareform(diss_matrix, force='tovector')
                print(f"The SNP list for this chromosome is: \n{snp_list}")
                print(f"With length of: {len(snp_list)}")
                result = 'OK'
        if i > 6:
            sys.exit('Invalid method.')
        if result == 'OK':
            if typed_input == 'hierarchy':
                hierarchy_clustering(diss_matrix_triangled, matrix_file_path)
            elif typed_input == 'dbscan':
                dbscan_clustering(diss_matrix, snp_list, matrix_file_path)
            elif typed_input == 'affinity_propagation':
                af_prop_clustering(corr_matrix, snp_list, matrix_file_path)
            elif typed_input == 'hdbscan':
                hdbscan_clustering(corr_matrix, snp_list, matrix_file_path)
            elif typed_input == 'k_means':
                kmeans_clustering(corr_matrix, snp_list, matrix_file_path)
            elif typed_input == 'k_medoids':
                kmedoids_clustering(corr_matrix, snp_list, matrix_file_path)
            elif typed_input == 'spectral':
                spectral_clustering(corr_matrix, snp_list, matrix_file_path)
            time_end = datetime.datetime.now()
            print(f"\nTime of finish - {time_end.isoformat(sep=' ', timespec='seconds')}. Time of executing - {time_end - time_start}.")
    except IndexError:
        sys.exit('Invalid method.')
    return met, matrix_file_path
awaiting_input()    
