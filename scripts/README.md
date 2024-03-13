# GF Detection Pipeline Requirements and Scripts
This section provides an overview of the programs, reference files, input files, and scripts required for the GF detection pipeline workflow. Each script has a brief description of its goal and expected output. We suggest creating a conda environment for running the programs and optimizing script runtime on the High-Performance Computing (HPC) cluster/computer for the best results.

## Programs 
```
Guppy6
Cutadapt (v4.6)
Minimap2 (v2.26)
Samtools (v1.19.2)
LongGF 
JAFFA (v2.3) 
FusionSeeker (v1.0.1) 
LongReadSum (v1.3.1)
Python (3.9)
```

## References
1. [Reference Genome: hg38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.15/)
2. [BED & GTF Comprehensive Gene Annotation from GENCODE](https://www.gencodegenes.org/human/) (best to use the most recent version)

## Input File Types
The pipeline is designed to take in a single input file, which can be in either `FAST5`, `FASTQ`, or `FASTA` format. If the input is in `FAST5` format, it is recommended to use Guppy6 (Available online: https://community.nanoporetech.com) to perform basecalling first, with the default basecalling read Q-score filtering threshold of 7, and then convert the FASTQ files for downstream analysis. If dealing with multiple `FASTQ` or `FASTA` files, it is best to merge the files using the `cat` command beforehand.  


## Overview 
![Rybacki_Fig1_pipeline_overview](https://github.com/WGLab/Gene-Fusion-Detection-Pipeline-LRS/assets/89222332/e4a3e393-126c-4b75-bff9-f0aebfbdce16)


## Scripts and Descriptions 
```
00_dir_setup.sh
01_redo_basecalling .sh
02_combine_cutadapt.sh
03_alignment_hg38.sh
04_longreadsum.sh
05_FusionSeeker.sh
05_JAFFA.sh
05_LongGF.sh
06_GF_criteria.py
```

### 0. Initialization of Directories
1. **Goal:** Generate the necessary folders to hold the downstream analysis results for each program.
2. **Expected Output:** 8 directories, following the numbering and names of the scripts
```
01_redo_basecalling
02_combine_cutadapt
03_alignment_hg38
04_longreadsum
05_JAFFA
05_LongGF
05_FusionSeeker
06_GF_criteria
```

### 1. Redo Basecalling with Guppy6
1. **Goal:** Redo the basecalling using a super high accuracy model from Guppy6.
2. **Expected Output:**







