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
    
    # for num_chr in range(13,23):
    num_chr = input('num_chr: ')
    with open(f"/mnt/wd/nsap/hap_assoc/chr{num_chr}_new.assoc.hap", 'r') as file:
        print(f"\r{num_chr} file...")
        pd_df1 = pd.read_csv(file, lineterminator='\n', delim_whitespace=True)
        print(pd_df1)
        pd_df = pd.concat([pd_df, pd_df1])
        # break
    print('DataFrame is merged.')
    
    ### sort dataframe by P-value
    pd_df = pd_df.sort_values(by=['P'])            
    ### creating and writing to the outfile
    ### OMNIBUS filtered
    omnibus_df = pd_df.where(pd_df['HAPLOTYPE'] == 'OMNIBUS').dropna(subset=['LOCUS'])
    with open(f"/mnt/wd/nsap/hap_assoc/csvs/chr{num_chr}_omnibus.csv", 'w') as omnibus_excel:
        omnibus_df.to_csv(omnibus_excel, line_terminator='\n', index=False)
        print('Created omnibus outfile...')
    non_omnibus_df = pd_df[pd_df.HAPLOTYPE != 'OMNIBUS']
    with open(f"/mnt/wd/nsap/hap_assoc/csvs/chr{num_chr}_nonOmnibus.csv", 'w') as non_omnibus_excel:
        non_omnibus_df.to_csv(non_omnibus_excel, line_terminator='\n', index=False)
        print('Created non omnibus outfile...')
    ### OMNIBUS not filtered
    # with open('old_sorted_chr1_.csv', 'w') as csv_file:
    #     pd_df.to_csv(csv_file, line_terminator='\n', index=False)
    print('Done!')
                    
dataframe_merger()
