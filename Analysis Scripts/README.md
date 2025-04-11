# Gene Fusion (GF) Detection and Analysis Pipeline

This pipeline provides a streamlined approach for identifying, combining, and filtering detected gene fusions using multiple long-read gene fusion detection tools and Jupyter Notebooks for analysis and result summarization. We suggest creating a conda environment for running the programs, `conda create -n GF_detection_LRS python=3.9`, and optimizing script runtime on the High-Performance Computing (HPC) cluster/computer for the best results.

### Programs 
```
Guppy6 (Dorado with recent updates)
Cutadapt (v4.6)
Minimap2 (v2.26)
Samtools (v1.19.2)
LongGF 
JAFFA (v2.3) 
FusionSeeker (v1.0.1) 
LongReadSum (v1.3.1)
Python (3.9)
```


## Overview

The figure below outlines the GF Detection and Analysis Pipeline with two distinct filtering strategies. Shared core components of the pipelines (center section, uncolored in the figure) are common to both strategies. Regardless of filtering strategy chosen, the pipeline begins by re-basecalling the long-read sequencing data using a super high-accuracy model (represented by the diamond symbol), Cutadapt if using the CHOP Cancer Fusion Panel approach, otherwise/then alignment to GRCh38 reference genome. From there, the analysis splits into two distinct filtering strategy paths. 

 
![Image](https://github.com/user-attachments/assets/d5d6c620-43ae-4e6e-987f-14facbff0a19)


### Input File Types
The pipeline is designed to take either `FAST5`, `FASTQ`, or `FASTA` file formats as input. If the input is in `FAST5` format, it is recommended to use Guppy6 (Available online: https://community.nanoporetech.com) to perform basecalling first, with the default basecalling read Q-score filtering threshold of 7, and then convert the FASTQ files for downstream analysis. If dealing with multiple `FASTQ` or `FASTA` files, it is best to merge the files using the `cat` command beforehand.  

### 1. Redo Basecalling with Guppy6 (Diamond Shape) 
1. **Goal:** Redo the basecalling using a super high accuracy model from Guppy6.
2. **Expected Output:** High-accuracy basecalled files, usually in the form of `FASTQ` or `FAST5` formats, distributed into corresponding "pass" and "fail" folders.



### 2. Combine and Cutadapt (Dependent on using the CHOP Cancer Fusion Panel Approach) 
1. **Goal:** Combine output `FASTQ` files from the "pass" folder of the re-basecalled samples and perform adapter trimming using the Cutadapt.
2. **Expected Output:** Final trimmed `FASTQ` file that contains trimmed reads from both adapter trimming processes (compressed in gzip format).



### 3. Alignment to the GRCh38 Reference Genome
1. **Goal:** Splice aware alignment to the hg38 reference using Minimap2.
2. **Expected Output:** `BAM` files containing alignment information, sorted by both read names and genomic position. 



### 4. Quality Check Analysis with LongReadSum 
1. **Goal:** Quality Check Analysis of the mapping to the hg38 reference with LongReadSum
2. **Expected Output:** `bam_summary.txt`, a summary file containing statistics of the BAM file and `QC_bam.html` HTML file with statistics of the BAM file.



### 5. Gene Fusion Detection: LongGF, JAFFAL, and FusionSeeker
1. **Goal:** Detect gene fusions; each program which can be run using their default settings. We recommended running **LongGF** with multiple parameters to optimize sensitivity and specificity. Below are the parameters used and corresponding output file name. 

| Output File Name          | Min Overlap | Bin Size | Min Map Len | 
|---------------------------|-------------|----------|-------------|
| longgf.100-25-100.txt     | 100         | 25       | 100         |
| longgf.100-10-100.txt     | 100         | 10       | 100         |
| longgf.75-25-75.txt       | 75          | 25       | 75          | 
| longgf.75-10-75.txt       | 75          | 10       | 75          | 
| longgf.50-25-50.txt       | 50          | 25       | 50          | 
| longgf.50-10-50.txt       | 50          | 10       | 50          | 
| longgf.25-25-25.txt       | 25          | 25       | 25          | 
| longgf.25-10-25.txt       | 25          | 10       | 25          | 
| longgf.10-25-10.txt       | 10          | 25       | 10          | 
| longgf.10-10-10.txt       | 10          | 10       | 10          | 

> **Note:**
> The pseudogene and Secondary_alignment parameters are set to 0 by default, meaning pseudogenes are not used from the GTF file, and secondary alignments are excluded during fusion detection.
> The min_sup_read parameter is set to 1 for all runs to detect all potential fusion events without any pre-filtering during fusion detection. Filtering is applied downstream during analysis and summarization.

2. **Expected Output:** Multiple fusion calls across the three programs. 

---



# GF Filtering 


### 1. **Combine Fusion Calls from All Programs**  
1. **Goal:** Summarize all fusion outputs from each fusion program individually and comprehensively composed of unique fusions based on fusion name, check the orientation of the fusion, and retain the entry with highest number of supporting reads if duplicated entries. *Notebook:* `01_Combine_Fusion_Calls.ipynb`  
2. **Expected Output:** A summarized fusion call file containing unique gene fusions (based on gene names and breakpoint positions) across all three programs with corresponding individual summary files from each program.


## 2. Filtering Strategies

#### A. CHOP Cancer Fusion Panel Filtering Pipeline (Top Section – Pink)
1. **Goal:** Filter gene fusions to focus on those involving genes listed in the CHOP Cancer Fusion Panel and remove likely artifactual fusions. Set a minimum read support threshold and recurrent percentage, and references provided. Notebook: `02_CHOP_Cancer_Panel_Filtering.ipynb`
2. **Expected Output:** A filtered list of clinically relevant fusions based on inclusion in the CHOP Cancer Panel and a count of remaining fusions after each filtering step.



#### B. Long-Read Whole Transcriptome Filtering Pipeline (Bottom Section – Orange)
1. **Goal:** Filter gene fusions to focus on those involving cancer-related genes, remove likely artifactual fusions, filter fusions present in parallel healthy tissue samples, and annotate remaining fusions based on their gene involvement in disease-specific fusions from the Mitelman database. Set a minimum read support threshold and recurrent percentage, and references provided. Notebook: `02_Long_Read_Transcriptome_Filtering.ipynb`
4. **Expected Output:** A filtered list of clinically relevant fusions based on inclusion in the cancer-related genes and a count of remaining fusions after each filtering step.




