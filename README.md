##SNP_clustering v0.5 

##10.06.2021

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

	Run the SNP_metadata_extractor.py script to create plink metadata.
	The script computes files called like "chr1_snplist.snplist" 
	and "chr1.ld". These files are created by SNP_metadata_extractor.py
	After all desired metadata is present run the snp_clustering_v0.5.py script. 
	Use the CLI to run the script with arguments:
	argv[1] being matrix file path
	argv[2] being snplist file path
	argv[3] being one of implemented clustering methods (dbscan, hdbscan, k-means, 
	k-medoids, spectral)
	affinity propagation and hierarchy.linkage are under tweaking and are advised not being used.
	
##Update

	tweaked (dbscan, hdbscan, k-means, k-medoids, spectral) parameters
	completed outfile creation in an optimised view

##WARNING

	This script is under development. Method parameters need to be tweaked for more precise results. 
	Please report any bugs to shushulba.rol@gmail.com
