# Mfuzz_RNAseq.R
A R script to perform clustering of gene expression time-series RNA-seq data with Mfuzz.

Required R libraries :
optparse,tools,Mfuzz,GenomicFeatures,DESeq,edgeR

Mfuzz webpage : http://mfuzz.sysbiolab.eu/
Mfuzz paper : http://w3.ualg.pt/%7Emfutschik/publications/bioinformation.pdf

Mfuzz_RNAseq.R take as input a set of RNA-seq count tables, one per sample, from HTSeq-count for example. 

Mfuzz_RNAseq.R performs a complete RNAseq data normalization and then uses Mfuzz package to perform a soft clustering of gene expression time-series data.

Normalization steps : From the input count tables, the Mfuzz_RNAseq.R script performs a library size normalization with DESeq method and then adjust these normalized data for gene length (normalized data / gene length). These normalization steps are carried out to make all the samples comparable, which is required by Mfuzz package.

Soft clustering steps : With these last normalized data (called RPKN data), the Mfuzz_RNAseq.R script performs a genes clustering analysis with Mfuzz package, generating clusters and associated genes lists.

This script has three principal inputs : 
- the argument "--folder" or "-f" which is the directory containing all the RNA-seq count tables (and only these files).
  Mfuzz_RNAseq.R will read and merge all these tables and will perform the normalization steps.
- the argument "--annotation" or "-a" is the path to an genes/transcripts annotation file (gff or gtf format), allowing 
  to calculate the genes length (sum of the exons length, overlap of exons is take into account). This lengths are used 
  during the data normalization by gene length.
- the argument "--time" or "-t" give the time value of each file by respecting the same order in the vector than the files in 
  the folder. This  is a list of type 'time1,time1,time1,time2,time2,time2,time3'. 
  If several files correspond to a same time (replicates), give the same time value and then the script performs the mean on 
  the normalized counts of all the samples of a same time to perform the soft clustering.

For a description of optional arguments, type : 
      /usr/bin/Rscript Mfuzz_RNAseq.R -h

Minimal command: 
      /usr/bin/Rscript Mfuzz_RNAseq.R -f count_files_folder -a annotation -t time

Complete command: 
      /usr/bin/Rscript Mfuzz_RNAseq.R -f count_files_folder -a annotation -b gene_name_attribute -t time -n nb_clusters 
      -m membership_cutoff -s min_std -e exclude_thres -r replacement_mode -o output_directory
