##SNP_clustering v0.7 

##18.08.2021

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
	mysql-connector

##Usage

	Run the SNP_metadata_extractor.py script to create plink metadata (SNP lists).
	The script computes files called like "chr1.snplist".
	After all desired metadata is present run the snp_clustering_v0.7_db.py script. 
	Use the CLI to run the script with arguments:
	argv[1] being matrix file path
	argv[2] being snplist file path
	argv[3] being one of implemented clustering methods (dbscan, hdbscan, k-means, spectral, all)
	k-medoids, affinity propagation and hierarchy.linkage are under tweaking and are not advised to use.
	Then to compare clustering results of working algorithms use cluster_MySQL_compare.py. This script
	will read data from DB (that has been created with the main clustering script in MySQL)
	and create files which are required as input for haplotype testing via Plink 1.07 
	(newer version of plink does not support haplotype testing, see [plink documentation](http://google.com))
	https://www.cog-genomics.org/plink2.
	cluster_MySQL_compare.py script has 2 functions which can be used separately but it's 
	recommended to run the 1st function (.hlist files creation) before 2nd (loop for executing 
	plink command that computes .frq.hap files for all existing infiles in the same directory).

	
##Update

	fully integrated with MySQL server for data storage
	fully working clusters comparison
	fully working infiles (for Plink 1.07 [command](https://zzz.bwh.harvard.edu/plink/haplo.shtml))
	tweaked (dbscan, hdbscan, k-means, k-medoids, spectral) parameters

##WARNING

	This script is under development. Method parameters need to be tweaked for more precise results. 
	Please report any bugs to shushulba.rol@gmail.com
