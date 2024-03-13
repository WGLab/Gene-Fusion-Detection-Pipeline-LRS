# Gene Fusion Detection Pipeline with Long-Read Sequencing 
The Gene Fusion (GF) Detection Pipeline is comprehensive bioinformatics workflow designed for GF detection from Oxford Nanopore long-read sequencing data. 

## Background 
Gene Fusions (GFs) occur when two distinct genes are fused together, forming a new hybrid gene. GFs typically result from genomic rearrangements including chromosomal inversions, translocations, duplications, and deletions. GFs often result in isoforms from the variable gene partners that attribute to the irregular protein coding genes, activation of oncogenes, or the inactivation of tumor suppressors, ultimately leading to cancer. GFs have been found in an increasing number of cancer related tumors and can be used as cancer biomarkers, as well as potential therapeutic targets. This pipeline leverages an ensemble of GF detection programs, namely, LongGF, JAFFA, and FusionSeeker, to identify GFs and incorporates filtering criteria, cross-referencing with established fusion databases, and manual validation. The GF Detection Pipeline was adapted from clinical cancer patient samples and works on long-read sequencing datasets, providing insights into known and potentially novel GFs, contributing to advancements in cancer diagnostics and precision medicine.

## Overview 
![Rybacki_Fig1_pipeline_overview](https://github.com/WGLab/Gene-Fusion-Detection-Pipeline-LRS/assets/89222332/e4a3e393-126c-4b75-bff9-f0aebfbdce16)

- **Black Arrows:** 29 de-identified clinical cancer patient samples previously processed by Illumina short-read sequencing library preparation and sequencing
- **Red Arrows:** Computational GF detection and analysis pipeline
- **Green Arrows:** Comparison to the results initially obtained by the DGD at CHOP
- **Blue Arrows:** Additional Potentially Novel GFs validated by PCR and Sanger Sequencing, followed by generation of fusion cDNA transcripts for functional validation in the Drosophila model for functional implications; potential incorporation into future GF analysis

## Program Requirements 
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

## Reference Requirements
```
[Reference Genome: hg38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.15/)
[Comprehensive gene annotation from GENCODE (BED & GTF)](https://www.gencodegenes.org/human/) (and/or the most recent version)
```

# Input File Types
The pipeline is designed to take in a single input file, which can be in either `FAST5`, `FASTQ`, or `FASTA` format. If the input is in `FAST5` format, it is recommended to use Guppy6 (Available online: https://community.nanoporetech.com) to perform basecalling first, with the default basecalling read Q-score filtering threshold of 7, and then convert the FASTQ files for downstream analysis. If dealing with multiple `FASTQ` or `FASTA` files, it is best to merge the files using the `cat` command beforehand.  

#### Note: Each example command may take some additional time ran as is, for best results, optimize with the HPC (High-Performance Computing) cluster/computer. 


## Script Descriptions 
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
2. **Example Command:** `./ 00_dir_setup.sh`
3. **Expected Output:** 8 directories, following the numbering and names of the scripts
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
2. **Example Command:** `./ 01_redo_basecalling .sh`
3. **Expected Output:** 






