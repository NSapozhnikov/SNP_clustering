import os

for num_chr in range(1, 23):
    os.system(f"plink --noweb --bfile merged --hap chr{num_chr}.hlist --hap-assoc --allow-no-sex --out chr{num_chr}") 
    os.remove(f"chr{num_chr}.log")
    os.remove(f"chr{num_chr}.mishap")
    os.remove(f"chr{num_chr}.nosex")
    print(f"{num_chr}/22...")
print('Done!')