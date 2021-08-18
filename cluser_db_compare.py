import pandas as pd
import numpy as np
import pyodbc
import os

cnxn = pyodbc.connect("Driver={SQL Server Native Client 11.0};"
                      "Server=DESKTOP-IHSSCUQ;"
                      "Database=SNP_clusters;"
                      "Trusted_Connection=yes;")
cursor = cnxn.cursor()
cursor.execute("""IF(NOT EXISTS(SELECT 1 FROM dbo.clusters_chr1))
                    BEGIN
                        RAISERROR('MyError',16,10);
                    END;""") 
cursor.commit()

# num_chr = input('Num_chr: ')
def preparing_infiles():
    for num_chr in range(1, 23):
        read_query = f"""SELECT SNP, dbscan, spectral, kmeans FROM dbo.clusters_chr{num_chr}
                            WHERE (dbscan NOT LIKE '-1');"""
        cursor.execute(read_query)
        clusters__list = []
        for row in cursor:
            clusters__list.append(row)
        print(clusters__list)
        k = 0
        fp = open(f"chr{num_chr}.hlist", 'w')
        for j in range(len(clusters__list)):    
            if j == 0:
                meta_list = []        
                meta_list.append(clusters__list[j][0])
                if (clusters__list[j][1], clusters__list[j][2], clusters__list[j][3]) == (clusters__list[j-1][1], clusters__list[j-1][2], clusters__list[j-1][3]):
                    meta_list.append(clusters__list[j][0])        
                else:       
                    for i in range(len(meta_list)-1):             
                        fp.write(f"** block{k} ")
                        for v in range(len(meta_list)):
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
                    if len(meta_list) > 1:
                        k+=1     
                        for i in range(len(meta_list)-1):             
                            fp.write(f"** block{k} ")
                            for v in range(len(meta_list)):
                                fp.write(' ')
                                fp.write(meta_list[v])
                            fp.write('\n')
                            meta_list.pop()
                        meta_list = []
                        meta_list.append(clusters__list[j][0])
        if len(meta_list) > 1:                
            k+=1
            for i in range(len(meta_list)-1):             
                fp.write(f"** block{k} ")
                for v in range(len(meta_list)):
                    fp.write(' ')
                    fp.write(meta_list[v])
                fp.write('\n')
                meta_list.pop()
        fp.close()
        
# preparing_infiles()        
def _hapfreq(): 
    print('Using plink...')
    i, j = 0, 0
    for num_chr in range(1, 23):
        try:
            os.system(f"plink --bfile merged --hap chr{num_chr}.hlist --hap-freq --out chr{num_chr}")
            os.remove(f"chr{num_chr}.log")
            os.remove(f"chr{num_chr}.mishap")
            os.remove(f"chr{num_chr}.nosex")
        except FileNotFoundError:
            i+=1       
        except FileExistsError:
            j+=1
        if i > 21:
            print('No infiles.')
        if j > 21:
            print("Infile already exists.")
    print('Done.')

# _hapfreq()    