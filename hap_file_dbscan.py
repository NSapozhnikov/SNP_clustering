import pandas as pd
import numpy as np
import datetime
import mysql.connector as mysql
from mysql.connector import ProgrammingError
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

def establish_engine():
    config = configparser.ConfigParser()
    config.read('credentials.ini')   
    engine = config['DEFAULT']['engine']
    return engine
  
connection = connectToBD_MySQL()
cursor = connection.cursor()
engine = create_engine(establish_engine())

num_chr = input('Num_chr: ')
### Update LOC column
try:
  upd_query = f"""UPDATE clusters_chr{num_chr}_new v
                  JOIN SNP_LOCATIONS_REP_CHR{num_chr} t ON v.SNP = t.SNP
                  SET v.LOC = t.LOC;"""
  with engine.begin() as conn:
      print('Joinning SNP locations...')
      conn.execute(upd_query)
      
  read_query = f"""SELECT SNP, dbscan, LOC
                   FROM clusters_chr{num_chr}_new
                   WHERE (dbscan NOT LIKE -1)
                   ORDER BY CAST(dbscan AS DECIMAL) ASC,
                            CAST(LOC AS DECIMAL) ASC;"""
  cursor.execute(read_query)
  clusters__list = cursor.fetchall()
  connection.commit()
  ### Create a dict with unique cluster labels
  counter = {}
  for i in range(len(clusters__list)):
    if clusters__list[i][2] != None:
      key = clusters__list[i][1]
      ## Append (SNP, LOC) tuple by unique keys
      if not key in counter:
          counter[key] = []
      counter[key].append([clusters__list[i][0], clusters__list[i][2]])
  ### Write to the outfile all values with unique keys with more than 1 element
  ### with hlist format for plink's 1.07 wildcard specification --hap flag
  with open(f"chr{num_chr}_new.hlist", 'w') as fp:
      clusters = 0
      for key in counter:
          if len(counter[key]) > 1:
              meta_list = counter[key]
              clusters+=1
              for j in range(len(meta_list)-1):
                  fp.write(f"** {num_chr}:")
                  snp = f"{meta_list[0][1]}-{meta_list[-1][1]}"
                  fp.write(snp)
                  for v in range(len(meta_list)):
                      fp.write(' ')
                      fp.write(meta_list[v][0])
                  fp.write('\n')
                  meta_list.pop()
      print(f"{clusters} clusters for chromosome {num_chr} has been written to the outfile.")
except ProgrammingError:
  print(f"There is no table for chromosome {num_chr} in the schema. Continuing...")
  
