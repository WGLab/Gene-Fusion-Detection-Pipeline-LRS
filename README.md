# Gene Fusion Detection Pipeline with Long-Read Sequencing 
The Gene Fusion (GF) Detection Pipeline is comprehensive bioinformatics workflow designed for GF detection from Oxford Nanopore long-read sequencing data. 

## Background 
Gene Fusions (GFs) occur when two distinct genes are fused together, forming a new hybrid gene. GFs typically result from genomic rearrangements including chromosomal inversions, translocations, duplications, and deletions. GFs often result in isoforms from the variable gene partners that attribute to the irregular protein coding genes, activation of oncogenes, or the inactivation of tumor suppressors, ultimately leading to cancer. GFs have been found in an increasing number of cancer related tumors and can be used as cancer biomarkers, as well as potential therapeutic targets. This pipeline leverages an ensemble of GF detection programs, namely, LongGF, JAFFAL, and FusionSeeker, to identify GFs and incorporates filtering criteria, cross-referencing with established fusion databases, and manual validation. The GF Detection Pipeline was adapted from clinical cancer patient samples and works on long-read sequencing datasets, providing insights into known and potentially novel GFs, contributing to advancements in cancer diagnostics and precision medicine.

## Overview 
![figure1_picture](https://github.com/user-attachments/assets/ae85e30b-168c-4e7b-8075-75a9d88fb12b)

- **Black Arrows:** 29 de-identified clinical cancer patient samples previously processed by Illumina short-read sequencing library preparation and sequencing
- **Red Arrows:** Computational GF detection and analysis pipeline
- **Green Arrows:** Comparison to the results initially obtained by the DGD at CHOP
- **Blue Arrows:** Additional Potentially Novel GFs validated by PCR and Sanger Sequencing

## Pipeline Setup and Script Details
For comprehensive information on setting up the GF detection pipeline workflow and understanding the scripts, navigate to the [scripts folder](https://github.com/WGLab/Gene-Fusion-Detection-Pipeline-LRS/tree/main/scripts). Further instructions and explanations are available in the dedicated README file.  

## Contact
If you have any questions, issues, or bugs, please post them on GitHub. 


