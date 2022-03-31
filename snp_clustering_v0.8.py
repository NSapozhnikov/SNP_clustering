"""snp_clustering v0.8
Please check README.txt"""

import numpy as np
import datetime
from sklearn.cluster import AffinityPropagation, DBSCAN, KMeans, SpectralClustering
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

def connectToBD_MySQL():
    config = configparser.ConfigParser()
    config.read('credentials.ini')    
    return mysql.connect(
            host="localhost",
            user=config['DEFAULT']['user'],
            password=config['DEFAULT']['password'],
            database='cl')
    

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
    db = DBSCAN(eps=0.7, min_samples=6, metric='precomputed').fit(diss_matrix)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
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
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine('mysql+mysqlconnector://nsap:iN5MrTmnP@localhost/cl')
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr} c
                    JOIN temporary_table u ON c.snp = u.snp
                    SET c.dbscan = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)

# def af_prop_clustering(corr_matrix, snp_list, matrix_file_path):
#     print('\nPerforming clustering with AffinityPropagation...')
#     clustering = AffinityPropagation(damping=0.99, max_iter=200,
#                                      convergence_iter=15, 
#                                      affinity='precomputed', random_state=0).fit(corr_matrix)
#     cluster_centers_indices = clustering.cluster_centers_indices_
#     labels = clustering.labels_
#     n_clusters_ = len(cluster_centers_indices)
#     print('Estimated number of clusters: %d' % n_clusters_)
#     out_file_path = matrix_file_path.rstrip('.ld')
#     with open(f"{out_file_path}_afprop.clustered", 'w') as fp:
#         print('Writting to the database...')
#         for i in range(n_clusters_):
#             cluster_indeces = np.where(labels == i)[0]
#             for j in cluster_indeces:
#                 #fp.write(str(snp_list[j]))
#                 #fp.write(' ')
#                 #fp.write(str(i))
#                 #fp.write('\n')

def hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with HDBSCAN...')
    clusterer = hdbscan.HDBSCAN(min_cluster_size=8, min_samples=6, metric='precomputed', cluster_selection_method='leaf').fit(diss_matrix)
    labels = clusterer.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
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
    print('Writting to the database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine('mysql+mysqlconnector://nsap:iN5MrTmnP@localhost/cl')
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr} c
                    JOIN temporary_table u ON c.snp = u.snp
                    SET c.hdbscan = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)

def kmeans_clustering(corr_matrix, snp_list, matrix_file_path, num_chr):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with K-means...')
    kmeans = KMeans(n_clusters=8, n_init=100).fit(corr_matrix)
    labels = kmeans.predict(corr_matrix)
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
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
    print('Writting to the database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine('mysql+mysqlconnector://nsap:iN5MrTmnP@localhost/cl')
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr} c
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
    clustering = SpectralClustering(n_clusters=8, affinity='precomputed').fit(corr_matrix)
    labels = clustering.labels_
    n_clusters_ = len(set(labels))
    print('Estimated number of clusters: %d' % n_clusters_)
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_spectral.clustered", 'w') as fp:
    # print('Writting to the database...')
    # for i in range(n_clusters_):
    #     cluster_indeces = np.where(labels == i)[0]
    #     for j in cluster_indeces:
    #         fp.write(str(snp_list[j]))
    #         fp.write(' ')
    #         fp.write(str(i))
    #         fp.write('\n')
    ### For output into DB
    print('Writting to the database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index':labels})
    engine = create_engine('mysql+mysqlconnector://nsap:iN5MrTmnP@localhost/cl')
    pd_df.to_sql(name='temporary_table', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE clusters_chr{num_chr} c
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
                        print('\nExtracting SNP list data from file...')
                        snp_list = np.loadtxt(snp_file, dtype=str)
                        cursor.execute(f"""SELECT * 
                                           FROM information_schema.tables
                                           WHERE table_schema = 'cl' 
                                           AND table_name = 'clusters_chr{num_chr}'
                                           LIMIT 1;""")
                        if cursor.fetchone() == None:
                            print('A new table will be created...')
                            cursor.execute(f"""CREATE TABLE clusters_chr{num_chr}(
                                                    SNP VARCHAR(32)
                                                    ,dbscan VARCHAR(32)
                                                    ,hdbscan VARCHAR(32)
                                                    ,spectral VARCHAR(32)
                                                    ,kmeans VARCHAR(32));""")
                        cursor.execute(f"""SELECT count(*) AS total FROM clusters_chr{num_chr};""") 
                        if cursor.fetchone()[0] == 0:
                            for snp in snp_list:
                                cursor.execute(f"""INSERT INTO cl.clusters_chr{num_chr} (SNP) VALUES ('{snp}')""")
                                connection.commit()
                except FileNotFoundError:
                    sys.exit('No such file.')
                try:
                    with open(matrix_file_path, 'r') as corr_file:
                        print('Extracting matrix data from file...')
                        corr = np.fromfile(corr_file, sep=' ')
                except FileNotFoundError:
                    sys.exit('No such file.')
                print('Reshaping an array...')
                corr = np.nan_to_num(corr)
                corr_matrix = corr.reshape(len(snp_list), len(snp_list))
                dissimilarity = 1 - abs(corr)
                diss_matrix = dissimilarity.reshape(len(snp_list), len(snp_list))
                # diss_matrix_triangled = squareform(diss_matrix, force='tovector')
                print(f"The SNP list for {str(num_chr)} chromosome is: \n{snp_list}")
                print(f"With length of: {len(snp_list)}")
                result = 'OK'
        if i > 7:
            sys.exit('Invalid method.')
        if result == 'OK':
            if typed_input == 'hierarchy':
                hierarchy_clustering(diss_matrix, matrix_file_path, num_chr)
            elif typed_input == 'dbscan':
                dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
            #elif typed_input == 'affinity_propagation':
            #    af_prop_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'hdbscan':
                hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'kmeans':
                kmeans_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'kmedoids':
                kmedoids_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'spectral':
                spectral_clustering(corr_matrix, snp_list, matrix_file_path, num_chr)
            elif typed_input == 'all':
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
