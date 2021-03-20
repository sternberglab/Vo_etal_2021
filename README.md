# Cointegrate analysis script

This script (`cointegrate_finder.py`) analyzes PacBio ccs reads to search for transposons, and classifies them according to whether they are from the originating plasmid, a simple genomic insertion, or a cointegrate. Full details are in the "SMRT data analysis" section of the manuscript, and comments are included in the script itself. 

## Requirements
This script uses bowtie2 and BLASTn to analyze and align reads. It also requires the python packages in the `requirements.txt` file. 

Much of the script depends on various input files (both samples being analyzed and associated fasta's for genomes, plasmids, and transposons) being accessible in an AWS S3 bucket. To adapt this for other data or purposes, point the file downloader to the location where those files are located in your setup. 

## Outputs
All outputs from the original analysis are included here for reference