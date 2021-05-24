# SNP_clustering v0.2
# Please check README.md

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
    
# """to do: implement data plotting"""

def hierarchy_clustering(diss_matrix_squared):
    print('Performing clustering with scipy.cluster.hierarchy...')
    print(f"Squared dissimilarity matrix is: {diss_matrix_squared}")
    print(f"With length of: {len(diss_matrix_squared)}")
    hierarchy = linkage(diss_matrix_squared, method='average')
    print('Output matrix is: ')
    print(hierarchy)
    fc = fcluster(hierarchy, 20, criterion='maxclust')
    plt.figure()
    dn = dendrogram(hierarchy, p=10, truncate_mode='level')
    plt.show()
# hierarchy_clustering()

def dbscan_clustering(diss_matrix, snp_list):
    print('Performing clustering with DBSCAN...')
    print('Using a dissimilarity matrix: ')
    print(diss_matrix)
    db = DBSCAN(eps=0.05, min_samples=2, metric='precomputed').fit(diss_matrix)
    # core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    # core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    for i in range(n_clusters_):
        cluster_indeces = np.where(labels == i)[0]
        print(f"Cluster {i} has {len(cluster_indeces)} values:")
        cluster_snp_list = []
        for j in cluster_indeces:
            cluster_snp_list.append(snp_list[j])
        if len(cluster_snp_list) < 3000:
            print(cluster_snp_list)
        else:
            print('This cluster contains the most SNPs.') 
# dbscan_clustering()

def af_prop_clustering(corr_matrix, snp_list):
    print('Performing clustering with AffinityPropagation...')
    print('Using a correlation matrix: ')
    print(corr_matrix)
    clustering = AffinityPropagation(damping=0.95, max_iter=4000,
                                     convergence_iter=400, random_state=None, 
                                     affinity='precomputed').fit(corr_matrix)
    cluster_centers_indices = clustering.cluster_centers_indices_
    labels = clustering.labels_
    n_clusters_ = len(cluster_centers_indices)
    print('Estimated number of clusters: %d' % n_clusters_)
    for i in range(n_clusters_):
        cluster_indeces = np.where(labels == i)[0]
        print(f"Cluster {i} has {len(cluster_indeces)} values:")
        cluster_snp_list = []
        for j in cluster_indeces:
            cluster_snp_list.append(snp_list[j])
        if len(cluster_snp_list) < 3000:
            print(cluster_snp_list)
        else:
            print('This cluster contains the most SNPs.')
# af_prop_clustering()

def hdbscan_clustering(corr_matrix): #doesn't build a minimum_spanning_tree with these parameters 
    print('Performing clustering with HDBSCAN...')
    print('Using a correlation matrix: ')
    print(corr_matrix)
    clusterer = hdbscan.HDBSCAN(metric='precomputed').fit(corr_matrix)
    # clusterer.minimum_spanning_tree_.plot(edge_cmap='viridis', edge_alpha=0.6, node_size=80, edge_linewidth=2)
    # clusterer.single_linkage_tree_.plot(cmap='viridis', colorbar=True)
    # clusterer.condensed_tree_.plot()
    # clusterer.condensed_tree_.plot(select_clusters=True, selection_palette=sns.color_palette())
# hdbscan_clustering()

def kmeans_clustering(corr_matrix, snp_list):
    print('Performing clustering with K-means...')
    kmeans = KMeans(n_clusters=5, n_init=100).fit(corr_matrix)
    print('Array of cluster affinity:')
    labels = kmeans.predict(corr_matrix)
    print(labels)
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    for i in range(n_clusters_):
        cluster_indeces = np.where(labels == i)[0]
        print(f"Cluster {i} has {len(cluster_indeces)} values:")
        cluster_snp_list = []
        for j in cluster_indeces:
            cluster_snp_list.append(snp_list[j])
        if len(cluster_snp_list) < 3000:
            print(cluster_snp_list)
        else:
            print('This cluster contains most SNPs.')
    trans = kmeans.transform(corr_matrix)
    print(trans)
# kmeans_clustering()
    
def kmedoids_clustering(corr_matrix, snp_list):
    print('Performing clustering with K-medoids...')
    kmedoids = KMedoids(n_clusters=8, metric='precomputed', max_iter=5).fit(corr_matrix)
    print('Array of cluster affinity:')
    labels = kmedoids.predict(corr_matrix)
    print(labels)
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    for i in range(n_clusters_):
        cluster_indeces = np.where(labels == i)[0]
        print(f"Cluster {i} has {len(cluster_indeces)} values:")
        cluster_snp_list = []
        for j in cluster_indeces:
            cluster_snp_list.append(snp_list[j])
        if len(cluster_snp_list) < 3000:
            print(cluster_snp_list)
        else:
            print('This cluster contains most SNPs.')
    trans = kmedoids.transform(corr_matrix)
    print(trans)    
# kmedoids_clustering()

def spectral_clustering(corr_matrix, snp_list):
    print('Performing clustering with SpectralClustering...')
    clustering = SpectralClustering(n_clusters=5, affinity='precomputed').fit(corr_matrix)
    labels = clustering.labels_
    print(labels)
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    for i in range(n_clusters_):
        cluster_indeces = np.where(labels == i)[0]
        print(f"Cluster {i} has {len(cluster_indeces)} values:")
        cluster_snp_list = []
        for j in cluster_indeces:
            cluster_snp_list.append(snp_list[j])
        if len(cluster_snp_list) < 3000:
            print(cluster_snp_list)
        else:
            print('This cluster contains most SNPs.')
# spectral_clustering()   

def awaiting_input():
    while True:
        try:
            typed_chr = int(input("Type in the chromosome number: "))
        except ValueError:
            print("Invalid input.")
            continue
        if typed_chr not in range(23):
            print("Must be in range [1; 22]")
            continue
        elif typed_chr == 0:
            print("Must be in range [1; 22]")
            continue
        else:
            break
    methods = ['hierarchy_clustering', 'dbscan_clustering', 'affinity_propagation_clustering',
               'hdbscan_clustering', 'k_means_clustering', 'k_medoids_clustering', 'spectral_clustering']
    print('Available methods: ')
    for mtd in methods:
        print(mtd)
    while True:
        result = 'NOT OK'
        i = 0
        typed_input = input('Type in the clustering method: ').lower().strip().replace(' ', '_').replace('-', '_')
        for met in methods:
            if typed_input not in met:
                i+=1                
            elif typed_input == '':
                i+=7         
            else:
                print(f"\n{met} algorithm will be used...")
                time_start = datetime.datetime.now()
                print (f"\nTime of start - {time_start.isoformat(sep=' ', timespec='seconds')}")    
                try:
                    with open(f"chr{typed_chr}_snplist.snplist", 'r') as snp_file:
                        print('\nExtracting SNP list data from file...')
                        snp_list = np.loadtxt(snp_file, dtype=str)
                except FileNotFoundError:
                    print('No such file.')
                    awaiting_input()
                try:
                    with open(f"chr{typed_chr}.ld", 'r') as corr_file:
                        print('Extracting matrix data from file...')
                        corr = np.fromfile(corr_file, sep=' ')
                except FileNotFoundError:
                    print('No such file.')
                    awaiting_input() 
                print('Reshaping an array...')
                corr_matrix = corr.reshape(len(snp_list), len(snp_list))
                dissimilarity = 1 - abs(corr)
                diss_matrix = dissimilarity.reshape(len(snp_list), len(snp_list))
                diss_matrix_squared = squareform(diss_matrix, force='tovector')
              # diss_matrix = pairwise_distances(corr_matrix)    
                print(f"The SNP list is: \n{snp_list}")
                print(f"With length of: {len(snp_list)}")
                result = 'OK' 
        if i > 6:
            print('Invalid method.')
        if result == 'OK':
            if typed_input in 'hierarchy_clustering':
                hierarchy_clustering(corr_matrix, snp_list)
            elif typed_input in 'dbscan_clustering':
                dbscan_clustering(diss_matrix, snp_list)
            elif typed_input in 'affinity_propagation_clustering':
                af_prop_clustering(corr_matrix, snp_list)
            elif typed_input in 'hdbscan_clustering':
                hdbscan_clustering(corr_matrix, snp_list)
            elif typed_input in 'k_means_clustering':
                kmeans_clustering(corr_matrix, snp_list)
            elif typed_input in 'k_medoids_clustering':
                kmedoids_clustering(corr_matrix, snp_list)
            elif typed_input in 'spectral_clustering':
                spectral_clustering(corr_matrix, snp_list)
            time_end = datetime.datetime.now()
            print(f"Time of finish - {time_end.isoformat(sep=' ', timespec='seconds')}. Time of executing - {time_end - time_start}.")
            break
       
awaiting_input()
