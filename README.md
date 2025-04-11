# Gene Fusion Detection Pipeline with Long-Read Sequencing 
The Gene Fusion (GF) Detection Pipeline is comprehensive bioinformatics workflow designed for GF detection from Oxford Nanopore long-read sequencing data analyzed with LongGF, JAFFAL, and FusionSeeker GF detection programs. 

## Background 
Gene Fusions (GFs) occur when two distinct genes are fused together, forming a new hybrid gene, that typically arrises from genomic rearrangements including chromosomal inversions, translocations, duplications, and deletions. GFs often result in isoforms that attribute to the irregular protein coding genes, activation of oncogenes, or the inactivation of tumor suppressors, ultimately leading to cancer. GFs have been found in an increasing number of cancer related tumors and can be used as cancer biomarkers, as well as therapeutic targets. 

This pipeline leverages an ensemble of GF detection programs (LongGF, JAFFAL, and FusionSeeker) to identify GFs and incorporates downstream filtering criteria, cross-referencing with established fusion databases, and computational validation. The GF Detection Pipeline was adapted from clinical cancer patient samples analyzed by the Children's Hospital of Philadelphia (CHOP) Cancer Fusion Panel and works on long-read sequenced datasets, providing insights into known and potentially novel GFs and contributing to advancements in cancer diagnostics and precision medicine.

## Overview of Sample Cohorts: 

| Sample Cohort | Description |
|--------|------------|
| **Cohort 1** | 29 short-read libraries previously sequenced with Illumina short-read sequencing, processed through the CHOP Cancer Fusion Panel, and later re-sequenced on the long-read ONT platform with Flongle flow cell. |
| **Cohort 1 Subset** | Same approach as Cohort 1, but with 24 of the 29 short-read libraries (Skipping Samples # 1, 3-6) and sequenced with ONT PromethION flow cell  |
| **Cohort 2** | 24 samples previously analyzed through the CHOP Cancer Fusion Panel and short-read sequencing, with total RNA prepared for ONT sequencing with PromethION flowcell. |

These two cohorts allow for a comparative analysis between short-read and long-read sequencing technologies in detecting GFs, providing insights into their respective strengths and limitations.


## Pipeline Setup and Script Details
For comprehensive information on using the GF detection pipeline analysis and filtering, navigate to the Analysis Scripts folder. Further instructions and explanations are available in the dedicated README file.  

## Contact
If you have any questions, issues, or bugs, please post them on GitHub. This work is part of the manuscript: _Combining Panel-Based and Whole Transcriptome Gene Fusion Detection Using Long-Read Sequencing_ (Manuscript in Preparation).
