import pandas as pd

def dataframe_merger():
    pd_df = pd.DataFrame({'LOCUS':[]
                         ,'HAPLOTYPE':[]
                         ,'F_A':[]
                         ,'F_U':[]
                         ,'CHISQ':[]
                         ,'DF':[]
                         ,'P':[]
                         ,'SNPS':[]})
    
    for num_chr in range(1,23):
        with open(f"chr{num_chr}.assoc.hap", 'r') as file:
            print(f"\r{num_chr} file...")
            for line in file:
                file_line = line.strip('\n')                                        ### remove the last symbol of the line
                file_line = file_line.split(sep=' ')                                ### split the line by space
                while '' in file_line: file_line.remove('')                         ### remove all remaining spaces
                if file_line[4] == 'CHISQ':                                         ### check for the header line
                    continue
                if file_line:            
                    pd_df1 = pd.DataFrame({'LOCUS': [file_line[0]]                  ### check if line is not empty -> append to dataframe
                                          ,'HAPLOTYPE':[file_line[1]]
                                          ,'F_A':[file_line[2]]
                                          ,'F_U':[file_line[3]]
                                          ,'CHISQ':[float(file_line[4])]
                                          ,'DF':[float(file_line[5])]
                                          ,'P':[float(file_line[6])]
                                          ,'SNPS':[file_line[7]]})
                    pd_df = pd_df.append(pd_df1, ignore_index=True)
                else:
                    print(num_chr, line)                                            ### print out number of file and 'empty' line 
    print('DataFrame is merged.')
    pd_df = pd_df.sort_values(by=['P'])                                             ### sort dataframe by P-value
    with open('out_merged_sorted.xlsx', 'w') as xlsx_file:
        pd_df.to_excel('out_merged_sorted.xlsx', index=False)                       ### creating and writing to the outfile
    print('Done!')
                    
