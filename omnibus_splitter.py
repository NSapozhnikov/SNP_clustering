import pandas as pd

file_path = 'J:/Visual Studio Projects/SNP_clustering/plink-1.07-dos/out_merged_sorted.xlsx'
with open(file_path, 'r') as excel_table:
    input_df = pd.read_excel(file_path)
    omnibus_df = input_df.where(input_df['HAPLOTYPE'] == 'OMNIBUS').dropna(subset=['LOCUS'])
    with open('out_merged_sorted_omnibus.xlsx', 'w') as omnibus_excel:
        omnibus_df.to_excel('out_merged_sorted_omnibus.xlsx', header=True, index=False)
        print('Created omnibus outfile...')
    non_omnibus_df = input_df[input_df.HAPLOTYPE != 'OMNIBUS']
    with open('out_merged_sorted_nonOmnibus.xlsx', 'w') as non_omnibus_excel:
        non_omnibus_df.to_excel('out_merged_sorted_nonOmnibus.xlsx', header=True, index=False)
        print('Created non omnibus outfile...')
print('Done!')