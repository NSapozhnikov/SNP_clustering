##SNP_clustering v0.4 

##29.05.2021

#SNP_clustering

	snp_clustering is a python script for clustering SNP correlation metadata.

##Requirements

		To use this script you need to have the following python packages:
	DateTime
	hdbscan 0.8.27
	NumPy 1.19.2
	matplotlib 3.3.2
	pandas 1.1.3
	scikit-learn 0.24.2
	scikit-learn-extra 0.2.0
	SciPy 1.6.2
	seaborn 0.11.1

##Usage

	Run the SNP_metadata_extractor.py script to create plink metadata. After all desired metadata is 
	present run the snp_clustering.py script. You will be prompted for entering the chromosome number 
	and clustering algorithm you wish to perform. The script computes files called like "chr1_snplist.snplist" 
	and "chr1.ld". These files are created by SNP_metadata_extractor.py
	
##Update

	Added out file creation
	Changed argument entry to bash

##WARNING

	This script is under development. Method parameters need to be tweaked for more precise results. 
	Please report any bugs to shushulba.rol@gmail.com
