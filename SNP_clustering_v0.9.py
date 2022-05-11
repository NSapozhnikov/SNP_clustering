"""snp_clustering v0.9
Please check README.txt"""

import numpy as np
import datetime
from sklearn.cluster import AffinityPropagation, DBSCAN, KMeans, SpectralClustering 
from sklearn.metrics.cluster import *
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from itertools import cycle
from sklearn.metrics import pairwise_distances
import hdbscan
from sklearn_extra.cluster import KMedoids
import pandas as pd
import sys
import argparse    
import mysql.connector as mysql
from getpass import getpass
import re
import configparser
from sqlalchemy import create_engine
from joblib import Memory


def connectToBD_MySQL():
    config = configparser.ConfigParser()
    config.read('credentials.ini')  
    return mysql.connect(
            host=config['DEFAULT']['host'],
            user=config['DEFAULT']['user'],
            password=config['DEFAULT']['password'],
            database='cl')

def establish_engine():
    config = configparser.ConfigParser()
    config.read('credentials.ini')   
    engine = config['DEFAULT']['engine']
    return engine

def hierarchy_clustering(diss_matrix, matrix_file_path):
    print('\nPerforming clustering with scipy.cluster.hierarchy...')
    hierarchy = linkage(diss_matrix, method='single')
    fcl_hierarchy = fcluster(hierarchy, 0.988, criterion='distance')
    print(fcl_hierarchy)
    print(fcl_hierarchy.max())

def dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with DBSCAN...')
    db = DBSCAN(eps=0.75, min_samples=6, metric='precomputed', n_jobs=-1).fit(diss_matrix)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    ### Scoring
    try:
        shil_score = silhouette_score(diss_matrix, labels, metric='precomputed') 
        print('Silhouette score: ', shil_score)
    except:
        print('Failed to compute Silhouette score.')
    try:
        CH = calinski_harabasz_score(diss_matrix, labels)
        print('Calinski-Harabasz score: ', CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_dbscan.clustered", 'w') as fp:
    # print('Writting to the outfile...)
    # for i in range(-1, n_clusters_):
    #     cluster_indices = np.where(labels == i)[0]
    #     for j in cluster_indices:
    #         fp.write(str(snp_list[j]))
    #         fp.write(' ')
    #         fp.write(str(i))
    #         fp.write('\n')
    ### For output into DB    
    print('\nWritting results to database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr}_new c
                    JOIN temporary_table u ON c.snp = u.snp
                    SET c.dbscan = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)

def af_prop_clustering(corr_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with AffinityPropagation...')
    clustering = AffinityPropagation(damping=0.99, max_iter=200,
                                     convergence_iter=15,
                                     affinity='precomputed', random_state=0).fit(corr_matrix)
    cluster_centers_indices = clustering.cluster_centers_indices_
    labels = clustering.labels_
    n_clusters_ = len(cluster_centers_indices)
    print('Estimated number of clusters: %d' % n_clusters_)


def hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
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
        shil_score = silhouette_score(diss_matrix, labels, metric='precomputed') 
        print('Silhouette score: ', shil_score)
    except:
        print('Failed to compute Silhouette score.')
    try:
        CH = calinski_harabasz_score(diss_matrix, labels)
        print('Calinski-Harabasz score: ', CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_hdbscan.clustered", 'w') as fp:  
    # print('Writting to the outfile...')
    # for i in range(-1, n_clusters_):
    #     cluster_indeces = np.where(labels == i)[0]
    #     for j in cluster_indeces:
    #         #fp.write(str(snp_list[j]))
    #         #fp.write(' ')
    #         #fp.write(str(i))
    #         #fp.write('\n')
    ### For output into DB 
    print('\nWritting results to database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr}_new c
                    JOIN temporary_table u ON c.snp = u.snp
                    SET c.hdbscan = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)
    mem.clear(warn=False)

def kmeans_clustering(corr_matrix, snp_list, matrix_file_path, num_chr):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with K-means...')
    kmeans = KMeans(n_clusters=1000, n_init=10, max_iter=100, n_jobs=-1).fit(corr_matrix)
    labels = kmeans.predict(corr_matrix)
    n_clusters_ = len(set(labels))
    print('Chosen number of clusters: %d' % n_clusters_)
    ### Scoring
    try:
        shil_score = silhouette_score(corr_matrix, labels) 
        print('Silhouette score: ', shil_score)
    except:
        print('Failed to compute Silhouette score.')
    try:
        CH = calinski_harabasz_score(corr_matrix, labels)
        print('Calinski-Harabasz score: ', CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_kmeans.clustered", 'w') as fp:
    # print('Writting to the outfile...')
    # for i in range(n_clusters_):
    #     cluster_indeces = np.where(labels == i)[0]
    #     for j in cluster_indeces:
    #         #fp.write(str(snp_list[j]))
    #         #fp.write(' ')
    #         #fp.write(str(i))
    #         #fp.write('\n')
    ### For output into DB
    print('\nWritting results to database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr}_new c
                    JOIN temporary_table u ON c.snp = u.snp
                    SET c.kmeans = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)
    
def kmedoids_clustering(corr_matrix, snp_list, matrix_file_path):
    print('\nPerforming clustering with K-medoids...')
    kmedoids = KMedoids(n_clusters=8, metric='correlation', method='pam', max_iter=1000).fit(corr_matrix)
    labels = kmedoids.predict(corr_matrix)
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_kmedoids.clustered", 'w') as fp:
    # print('Writting to outfile...')
    # for i in range(n_clusters_):
    #     cluster_indeces = np.where(labels == i)[0]
    #     for j in cluster_indeces:
    #         fp.write(str(snp_list[j]))
    #         fp.write(' ')
    #         fp.write(str(i))
    #         fp.write('\n')
                
def spectral_clustering(corr_matrix, snp_list, matrix_file_path, num_chr):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with SpectralClustering...')
    clustering = SpectralClustering(n_clusters=1000, affinity='precomputed', n_jobs=-1).fit(corr_matrix)
    labels = clustering.labels_
    n_clusters_ = len(set(labels))
    print('Chosen number of clusters: %d' % n_clusters_)
    ### Scoring
    try:
        shil_score = silhouette_score(corr_matrix, labels)
        print('Silhouette score: ', shil_score)
    except:
        print('Failed to compute Silhouette score.')
    try:
        CH = calinski_harabasz_score(corr_matrix, labels)
        print('Calinski-Harabasz score: ', CH)
    except:
        print('Failed to compute Calinski-Harabasz score.')
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_spectral.clustered", 'w') as fp:
    # print('Writting to the outfile...')
    # for i in range(n_clusters_):
    #     cluster_indeces = np.where(labels == i)[0]
    #     for j in cluster_indeces:
    #         fp.write(str(snp_list[j]))
    #         fp.write(' ')
    #         fp.write(str(i))
    #         fp.write('\n')
    ### For output into DB
    print('\nWritting results to database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr}_new c
                    JOIN temporary_table u ON c.snp = u.snp
                    SET c.spectral = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)

def awaiting_input():
    try:
        matrix_file_path = str(sys.argv[1])
        f = open(matrix_file_path, 'r')
        f.close()
        num_chr = re.findall("\\D*(\\d*)", matrix_file_path)[0]       
    except (IndexError, FileNotFoundError, ValueError ):
        sys.exit('Invalid matrix file path.')
    try:
        snp_file_path = str(sys.argv[2])
        f = open(snp_file_path, 'r')
        f.close()
        num_chr1 = re.findall("\\D*(\\d*)", snp_file_path)[0]
        if num_chr1 != num_chr:
            sys.exit('Chromosome missmatch in file names.')
    except (IndexError, FileNotFoundError, ValueError ):
        sys.exit('Invalid SNP list file path.')    
    try:
        methods = ['hierarchy', 'dbscan', 'affinity_propagation',
               'hdbscan', 'kmeans', 'kmedoids', 'spectral', 'all']
        result = 'NOT OK'
        i = 0
        typed_input = sys.argv[3].lower().replace('-', '_')
        for met in methods:
            if typed_input != met:
                i+=1                
            elif typed_input == '':
                i+=8         
            else:
                print(f"\n{met} algorithm will be used...")
                time_start = datetime.datetime.now()
                print (f"\nTime of start - {time_start.isoformat(sep=' ', timespec='seconds')}") 
                connection = connectToBD_MySQL()
                cursor = connection.cursor()
                try:
                    with open (snp_file_path, 'r') as snp_file: 
                        engine = create_engine(establish_engine())
                        print('\nExtracting SNP list data from file...')
                        snp_list = np.loadtxt(snp_file, dtype=str)        
                        print(f"The SNP list for {str(num_chr)} chromosome is: \n{snp_list}")
                        print(f"With length of: {len(snp_list)}")
                        cursor.execute(f"""SELECT * 
                                           FROM information_schema.tables
                                           WHERE table_schema = 'cl' 
                                           AND table_name = 'clusters_chr{num_chr}_new'
                                           LIMIT 1;""")
                        if cursor.fetchone() == None:
                            print(f"No table in the schema was found for chromosome {num_chr}.", 
                                     'A new table will be created...', sep='\n')
                            cursor.execute(f"""CREATE TABLE clusters_chr{num_chr}_new(
                                                    SNP VARCHAR(32)
                                                    ,LOC VARCHAR(32)
                                                    ,dbscan VARCHAR(32)
                                                    ,hdbscan VARCHAR(32)
                                                    ,spectral VARCHAR(32)
                                                    ,kmeans VARCHAR(32));""")
                            connection.commit()
                        cursor.execute(f"""SELECT count(*) AS total FROM clusters_chr{num_chr}_new;""")
                        if cursor.fetchone()[0] == 0:
                            baseSQL = f"INSERT INTO cl.clusters_chr{num_chr}_new (SNP) VALUES"
                            for snp in snp_list:
                                baseSQL += "('" + snp + "'),"
                            baseSQL = baseSQL.rstrip(',')
                            with engine.begin() as conn:
                                conn.execute(baseSQL)
                except FileNotFoundError:
                    sys.exit('No such file.')
                try:
                    with open(matrix_file_path, 'r') as corr_file:
                        print('Extracting matrix data from file...')
                        corr_matrix = np.fromfile(corr_file, sep=' ')
                except FileNotFoundError:
                    sys.exit('No such file.')
                np.nan_to_num(corr_matrix, copy=False)
                # diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
                # diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
                # np.fill_diagonal(diss_matrix, 0)
                # corr_matrix = corr_matrix.reshape(len(snp_list), len(snp_list))
                # np.fill_diagonal(corr_matrix, 1)                
                result = 'OK'
        if i > 7:
            sys.exit('Invalid method.')
        if result == 'OK':
            if typed_input == 'hierarchy':
                print('Reshaping an array...')
                diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
                corr_matrix = None
                diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
                np.fill_diagonal(diss_matrix, 0)
                hierarchy_clustering(diss_matrix, matrix_file_path, num_chr)
            elif typed_input == 'dbscan':
                print('Reshaping an array...')
                diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
                corr_matrix = None
                diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
                np.fill_diagonal(diss_matrix, 0)
                dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
            #elif typed_input == 'affinity_propagation':
            #    af_prop_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'hdbscan':
                print('Reshaping an array...')
                diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
                corr_matrix = None
                diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
                np.fill_diagonal(diss_matrix, 0)              
                hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'kmeans':
                print('Reshaping an array...')
                corr_matrix = corr_matrix.reshape(len(snp_list), len(snp_list))
                np.fill_diagonal(corr_matrix, 1)
                kmeans_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            # elif typed_input == 'kmedoids':
            #     print('Reshaping an array...')
            #     corr_matrix = corr_matrix.reshape(len(snp_list), len(snp_list))
            #     np.fill_diagonal(corr_matrix, 1)              
            #     kmedoids_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'spectral':
                print('Reshaping an array...')
                corr_matrix = corr_matrix.reshape(len(snp_list), len(snp_list))
                np.fill_diagonal(corr_matrix, 1)              
                spectral_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'all':
            ### Keeps both simmilarity and dissimilarity matrices in the RAM,
            ### which is twice RAM needed
                print('Warning: may have a high impact on RAM!')
                print('Reshaping arrays...')
                diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
                diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
                np.fill_diagonal(diss_matrix, 0)
                corr_matrix = corr_matrix.reshape(len(snp_list), len(snp_list))
                np.fill_diagonal(corr_matrix, 1)
                dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
                kmeans_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
                spectral_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)  
                hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
            time_end = datetime.datetime.now()
            print(f"""\nTime of finish - {time_end.isoformat(sep=' ', timespec='seconds')}. 
                      \nTime of executing - {time_end - time_start}.""")
    except IndexError:
        sys.exit('Invalid method.')
    return met, matrix_file_path
awaiting_input()   
