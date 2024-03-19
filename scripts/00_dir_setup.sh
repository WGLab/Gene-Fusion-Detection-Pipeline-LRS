#!/bin/bash

# Store this script in the main folder where you would like to generate the folder for downstream analysis 
# All folders will be generated, however depending on your input file type you may skip some steps, as described below:


# If your input is long-read sequenced reads from Nanopore sequencing:

# 1. Re-do Basecalling with Guppy6
mkdir -p "01_redo_basecalling"

# 2. Combine the Passed reads & Illumina adapter trimming for analysis 
mkdir -p "02_combine_cutadapt"


# If you have long reads in fastq format:

# 3. Align long reads to the reference genome of interest 
mkdir -p "03_alignment_hg38"


# LongReadSum Quality Check Analysis
mkdir -p "04_longreadsum"


# 5. Gene Fusion Detection 
mkdir -p "05_JAFFA"
mkdir -p "05_LongGF"
mkdir -p "05_FusionSeeker" 

# 6. Gene Fusion Filtering Criteria and Summary Files
mkdir -p "06_GF_criteria"


echo "Finished generating directories for downstream analysis" 


