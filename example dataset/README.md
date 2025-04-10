# Gene Fusion Detection and Analysis Pipeline

This pipeline provides a streamlined approach for identifying, combining, and filtering detected gene fusions using multiple long-read gene fusion detection tools and Jupyter Notebooks for analysis and result summarization.

## Overview

The pipeline consists of **three core fusion detection programs**, which can be run using their **default settings** for standard usage:

- [LongGF](https://github.com/WGLab/LongGF)  
- [JAFFAL](https://github.com/Oshlack/JAFFA/wiki)  
- [FusionSeeker](https://github.com/Maggi-Chen/FusionSeeker)

**We recommended running LongGF with multiple parameters to optimize sensitivity and specificity.** Below are the parameters used and corresponding output file name. 

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


---


## 1. **Combine All Program Fusion Calls**  

**Description:**
Summarize all fusion outputs from each fusion program individually and comprehensively composed of unique fusions based on fusion name, check the orientation of the fusion, and retain the entry with highest number of supporting reads if duplicated entries. 
> *Notebook:* `01_Combine_Fusion_Calls.ipynb`  



---



## 2. **Gene Fusion Detection and Analysis Pipelines**

The figure below outlines two distinct filtering strategies that can be applied following gene fusion detection.
![Image](https://github.com/user-attachments/assets/d5d6c620-43ae-4e6e-987f-14facbff0a19)

Both pipelines begin by re-basecalling the long-read sequencing data using a super high-accuracy model (represented by the diamond symbol). From there, the analysis splits into two distinct filtering strategy paths. 
Shared core components of the pipelines (center section, uncolored in the figure) are common to both strategies. These steps are indicated by black arrows and are applied to all samples regardless of the filtering path chosen.


### Filtering Strategies:

1. **CHOP Cancer Fusion Panel Filtering Pipeline** (top section, pink)
   - *Notebook:* `02_CHOP_Cancer_Panel_Filtering.ipynb`  
   - Used for targeted fusion detection in cancer-specific genes

2. **Long-Read Whole Transcriptome Filtering Pipeline** (bottom section, orange)
   - *Notebook:* `02_Long_Read_Transcriptome_Filtering.ipynb`  
   - Used for unbiased, comprehensive fusion discovery across the whole transcriptome 

You can choose either strategy based on your filtering preferences:
- **Option 1** (CHOP Cancer Fusion Panel Genes, pink): Focuses on fusions involving genes targeted in the CHOP Cancer Fusion Panel, using a minimum supporting read threshold of 2 or 10.
- **Option 2** (Whole Transcriptome, orange): Focuses on fusions involving cancer-related genes, filters out fusions present in parallel healthy samples, and annotates remaining fusions based on their gene involvement in disease-specific fusions from the Mitelman database. 

---






