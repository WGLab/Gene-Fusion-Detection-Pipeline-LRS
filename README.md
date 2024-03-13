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


## Contact
If you have any questions/issues/bugs, please post them on GitHub. 


