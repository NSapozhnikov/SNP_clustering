import pandas as pd
import mne

path = 'J:/Visual Studio Projects/SNP_clustering/'
omnibus_pd = pd.read_excel(path + 'out_merged_sorted_omnibus.xlsx')
nonOmnibus_pd = pd.read_excel(path + 'out_merged_sorted_nonOmnibus.xlsx')       ### Reading files

# fdr_omnibus = mne.stats.fdr_correction(omnibus_pd.P)                                
# fdr_nonOmnibus = mne.stats.fdr_correction(nonOmnibus_pd.P)                      ### Evaluating FDR

# omnibus_pd['FDR'] = fdr_omnibus[1]                                          
# nonOmnibus_pd['FDR'] = fdr_nonOmnibus[1]

bonferroni_omnibus = mne.stats.bonferroni_correction(omnibus_pd.P)
bonferroni_nonOmnibus = mne.stats.bonferroni_correction(nonOmnibus_pd.P)

omnibus_pd['bonferroni'] = bonferroni_omnibus[1]                                          
nonOmnibus_pd['bonferroni'] = bonferroni_nonOmnibus[1]                                        ### Adding up to the main dataframes   

with open(path + 'omnibus1.xlsx', 'w') as omni_fp:
    omnibus_pd.to_excel(path + 'omnibus1.xlsx', header=True, index=False)
    print('Writing omnibus.xlsx')
with open(path + 'nonOmnibus1.xlsx', 'w') as nonomni_fp:
    nonOmnibus_pd.to_excel(path + 'nonOmnibus1.xlsx', header=True, index=False)
    print('Writing nonOmnibus.xlsx...')                                         ### Writing outfiles
print('Done!')