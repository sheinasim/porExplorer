# porExplorer
Visualize your raw Nanopore or PacBio sequences from fastq

environment: python 2
requires: numpy and matplotlib
usage: porExplorer.py <fastq> <outfile prefix> <genome size in Mb (write 100 Mb genome as 100)> <desired coverage (write 40x as 40)>

outputs: 
1) accumulation plot showing coverage and minimum sequence size to achieve desired coverage
2) a histogram of raw read lengths
3) a histogram of read lengths as a proportion of the data
