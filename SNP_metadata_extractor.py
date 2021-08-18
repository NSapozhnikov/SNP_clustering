#!/usr/bin/env python
# coding: utf-8

# In[19]:


import os

for i in range(22):
    os.system(f"./plink --bfile merged --chr {i+1} --write-snplist --out chr{i+1}")
    os.system(f"rm chr{i+1}.log chr{i+1}.nosex")
    #chr{i+1}.ld
    print(f"{i+1} element has been computed...")
print('Done')


# In[ ]:




