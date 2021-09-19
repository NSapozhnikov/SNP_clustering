import pandas as pd
import numpy as np
import pyodbc
import os
import mysql.connector as mysql
from mysql.connector import ProgrammingError
import configparser

# cnxn = pyodbc.connect("Driver={SQL Server Native Client 11.0};"
#                       "Server=DESKTOP-IHSSCUQ;"
#                       "Database=SNP_clusters;"
#                       "Trusted_Connection=yes;")
# cursor = cnxn.cursor()
# cursor.execute("""IF(NOT EXISTS(SELECT 1 FROM dbo.clusters_chr1))
#                     BEGIN
#                         RAISERROR('MyError',16,10);
#                     END;""") 
# cursor.commit()
config = configparser.ConfigParser()
config.read('credentials.ini')    
connection = mysql.connect(
                host="localhost",
                user=config['DEFAULT']['user'],
                password=config['DEFAULT']['password'],
                database='SNP_Clusters')
cursor = connection.cursor()

# def create_DataFrame():
    # with open(f"J:/Visual Studio Projects/SNP_clustering/plink-1.07-dos/merged.bim", 'r') as file:
    #     pd_df = pd.DataFrame({'CHR': []                 
    #                           ,'SNP':[]
    #                           ,'LOCATION':[]})                                   ### This block is for the 1st use to read the .bim file 
    #     for line in file:                                                        ### with SNP locations. It will create .pkl file to read DataFrame from
    #             file_line = line.strip('\n')                                     ### remove the last symbol of the line
    #             file_line = file_line.split(sep='\t')                            ### split the line by space
    #             while '' in file_line: file_line.remove('')                      ### remove all remaining spaces
    #             if file_line: 
    #                 print(file_line)
    #                 pd_df1 = pd.DataFrame({'CHR': [file_line[0]]                 ### check if line is not empty -> append to dataframe
    #                                       ,'SNP':[file_line[1]]
    #                                       ,'LOCATION':[file_line[3]]})
    #                 pd_df = pd_df.append(pd_df1, ignore_index=True)
    # with open('pd_df_store.pkl', 'w') as f:
    #     pd_df.to_pickle('pd_df_store.pkl')
# create_DataFrame()

def preparing_infiles():
    pd_df = pd.read_pickle('pd_df_store.pkl')                                  ### This DataFrame is created with the Create_DataFrame() function 
    for num_chr in range(1, 23):
        try:
            read_query = f"""SELECT SNP, dbscan, spectral, kmeans FROM clusters_chr{num_chr}
                                WHERE (dbscan NOT LIKE '-1');"""
            cursor.execute(read_query)
            clusters__list = []
            for row in cursor:                                                 ### Reading all rows from  the table in DB to the list
                clusters__list.append(row)
            print(f"{num_chr}/22...")
            fp = open(f"chr{num_chr}.hlist", 'w')                              ### For the 1st element
            for j in range(len(clusters__list)):                               ### Comparing values in columns of the same row (different clustering algorithms)
                if j == 0:
                    meta_list = []        
                    meta_list.append(clusters__list[j][0])
                    if (clusters__list[j][1], clusters__list[j][2], clusters__list[j][3]) == (clusters__list[j-1][1], clusters__list[j-1][2], clusters__list[j-1][3]):
                        meta_list.append(clusters__list[j][0])                 ### adding the 1st SNP
                    else:                                                      
                        for i in range(len(meta_list)-1):                             
                            fp.write(f"** {num_chr}:")
                            snp = f"{pd_df.loc[pd_df['SNP'] == meta_list[0],'LOCATION'].item()}-{pd_df.loc[pd_df['SNP'] == meta_list[-1],'LOCATION'].item()}"
                            fp.write(snp)                                 
                            for v in range(len(meta_list)):                    ### adding all other SNPs ** chr:startSNP-stopSNP rs...                 
                                fp.write(' ')
                                fp.write(meta_list[v])
                            fp.write('\n')
                            meta_list.pop()
                        meta_list = []
                        meta_list.append(clusters__list[j][0])
                else:
                    if (clusters__list[j][1], clusters__list[j][2], clusters__list[j][3]) == (clusters__list[j-1][1], clusters__list[j-1][2], clusters__list[j-1][3]):
                        meta_list.append(clusters__list[j][0])        
                    else:
                        if len(meta_list) > 1:                                 ### same for other elements  
                            for i in range(len(meta_list)-1):             
                                fp.write(f"** {num_chr}:")
                                snp = f"{pd_df.loc[pd_df['SNP'] == meta_list[0],'LOCATION'].item()}-{pd_df.loc[pd_df['SNP'] == meta_list[-1],'LOCATION'].item()}"
                                fp.write(snp)
                                for v in range(len(meta_list)):                                    
                                    fp.write(' ')
                                    fp.write(meta_list[v])
                                fp.write('\n')
                                meta_list.pop()
                            meta_list = []
                            meta_list.append(clusters__list[j][0])
            if len(meta_list) > 1:                
                for i in range(len(meta_list)-1):             
                    fp.write(f"** {num_chr}:")
                    snp = f"{pd_df.loc[pd_df['SNP'] == meta_list[0],'LOCATION'].item()}-{pd_df.loc[pd_df['SNP'] == meta_list[-1],'LOCATION'].item()}"
                    fp.write(snp)
                    for v in range(len(meta_list)):                        
                        fp.write(' ')
                        fp.write(meta_list[v])
                    fp.write('\n')
                    meta_list.pop()
            fp.close()
        except ProgrammingError:
            continue
        finally:
            print('Done!')
preparing_infiles()   
     
# def _hapfreq(): 
#     print('Using plink...')
#     i, j = 0, 0
#     for num_chr in range(1, 23):
#         try:
#             os.system(f"plink --bfile merged --hap chr{num_chr}.hlist --hap-freq --out chr{num_chr}")
#             os.remove(f"chr{num_chr}.log")
#             os.remove(f"chr{num_chr}.mishap")                                                         ###
#             os.remove(f"chr{num_chr}.nosex")                                                          ###
#         except FileNotFoundError:                                                                     ###
#             i+=1                                                                                      ###                
#         except FileExistsError:                                                                       ### Depricated. Use assoc_hap.py
#             j+=1                                                                                      ###    
#         if i > 21:                                                                                    ###    
#             print('No infiles.')                                                                      ###            
#         if j > 21:                                                                                    ###
#             print("Infile already exists.")
#     print('Done.')

# # _hapfreq()    
