# GF Detection Pipeline Requirements and Scripts
This section provides an overview of the programs, reference files, input files, and scripts required for the GF detection pipeline workflow. Each script has a brief description of its goal and expected output. We suggest creating a conda environment for running the programs, `conda create -n GF_detection_LRS python=3.9`, and optimizing script runtime on the High-Performance Computing (HPC) cluster/computer for the best results.

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
3. CHOP Fusion Panel Genes
4. Mitelman Gene Fusion Database
5. COSMIC Gene Fusion Database
6. ChimerDB Gene Fusion Database


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

## 0. Initialization of Directories `00_dir_setup.sh`
1. **Goal:** Generate the necessary folders to hold the downstream analysis results for each program.
2. **Example Command:** `bash 00_dir_setup.sh`
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

## 1. Redo Basecalling with Guppy6 `01_redo_basecalling.sh`
1. **Goal:** Redo the basecalling using a super high accuracy model from Guppy6.
2. **Example Command:** `bash 01_redo_basecalling.sh sample_name input_dir output_dir basecaller_input config_file num_callers gpu_runners_per_device cpu_threads_per_caller device_cuda`
- `sample_name`: Enter the name for the sample.
- `input_dir`: Specify the path to the input directory.
- `output_dir`: Specify the full path to 01_redo_basecalling output. 
- `basecaller_input`: Specify the path to the basecaller model. 
- `config_file`: Specify the path to the configuration file name.
- Additional HPC Machine Configuration (Adjust based on your HPC machine specifications):
   - `num_callers`, `gpu_runners_per_device`, `cpu_threads_per_caller`, and `device_cuda` 

3. **Expected Output:** High-accuracy basecalled files, usually in the form of `FASTQ` or `FAST5` formats, distributed into corresponding "pass" and "fail" folders.






## 2. Combine and Cutadapt `02_combine_cutadapt.sh`
1. **Goal:** Combine output `FASTQ` files from the "pass" folder of the re-basecalled samples and perform adapter trimming using the Cutadapt.
2. **Example Command:** `bash 02_combine_cutadapt.sh sample_name pass_input_dir output_dir adapter3 adapter5 rev_adapter3 rev_adapter5 error_rate_trim1 error_rate_trim2 num_occurances_trim1 num_occurances_trim2`
- `sample_name`: Enter the desired name for the sample.
- `pass_input_dir`: Specify the full path to the "pass" folder of the re-basecalled samples.
- `output_dir`: Specify the full path to 02_combine_cutadapt output. 
- `adapter3`, `adapter5`, `rev_adapter3`, and `rev_adapter5`: Forward and Reverse adapter sequences for trimming.
- `error_rate_trim1` and `error_rate_trim2`: Cutadapt error rate denoted for both trimming steps. 
- `num_occurances_trim1` and `num_occurances_trim2`: Reads with fewer than the specified occurrences of a specified adapter were trimmed, denoted for both trimming steps. 

3. **Expected Output:** See below.

| Output File Names      | Description |
| -------------- | ----------- |
| `sample_name_pass_reads.fastq`            | A consolidated `FASTQ` file containing the passed reads after re-basecalling, meeting the specified quality criteria.       |
| `sample_name_pass_reads.fastq.gz`         | Similar to the above, but the file is compressed in gzip format (gz).        |
| `sample4_pass_reads_trimmed1.fastq`       | Trimmed `FASTQ` file of reads that were trimmed after the initial trimming process; adapter removal around the specified region of interest. | 
| `sample4_pass_reads_untrimmed1.fastq`     | Untrimmed `FASTQ` file containing reads that were not trimmed during the initial process; adapter removal around the specified region of interest. | 
| `cutadapt_report_linked_adapters.txt`     | A report generated by Cutadapt, detailing information about the adapter trimming process, particularly for linked adapters. | 
| `sample4_pass_reads_trimmed2.fastq`       | Trimmed `FASTQ` file of reads that were trimmed after the initial trimming process; adapter removal based on the specified individual adapter sequences.| 
| `sample4_pass_reads_untrimmed2.fastq`     | Untrimmed `FASTQ` file containing reads that were not trimmed during the initial process; adapter removal based on the specified individual adapter sequences.| 
| `cutadapt_report_single_adapters.txt`       | A report generated by Cutadapt, detailing information about the trimming process for single adapters. |
| `sample4_pass_reads_final_trimmed.fastq`    | Final trimmed `FASTQ` file; contains trimmed reads from both adapter trimming processes. |
| `sample4_pass_reads_final_trimmed.fastq.gz` | Similar to the above, but the file is compressed in gzip format (gz). | 






## 3. Alignment to the hg38 Reference `03_alignment_hg38.sh`
1. **Goal:** Splice aware alignment to the hg38 reference using Minimap2.
2. **Example Command:** `bash 03_alignment_hg38.sh sample_name input_fastq output_dir REF_GENOME REF_SPLICE`
- `sample_name`: Enter the name for the sample.
- `input_fastq`: Specify the full path to the final trimmed FASTQ file.
- `output_dir`: Specify the full path to 03_alignment_hg38 output. 
- `REF_GENOME`: Specify the full path to the hg38 reference genome file (in FASTA format).
- `REF_SPLICE`: Specify the full path to the BED file containing splice site annotations for hg38.

3. **Expected Output:** See below.

| Output File Names      | Description |
| -------------- | ----------- |
| `<sample_name>-guppy6.bam`                        | BAM file containing alignment information.        |
| `<sample_name>-guppy6.sam`                        | SAM file containing alignment information.        |
| `<sample_name>-guppy6.sorted_name.bam`            | BAM file sorted by read names.                    |
| `<sample_name>-guppy6.sorted_position.bam`        | BAM file sorted by genomic position.              |
| `<sample_name>-guppy6.sorted_position.bam.bai`    | Index file for the sorted BAM file by position.   |






## 4. Quality Check Analysis with LongReadSum  `04_longreadsum.sh`
1. **Goal:** Quality Check Analysis of the mapping to the hg38 reference with LongReadSum
2. **Example Command:** `bash 04_longreadsum.sh sample_name output_dir input_BAM long_read_sum_path`
- `sample_name`: Enter the name for the sample.
- `output_dir`:  Specify the full path to 04_longreadsum output. 
- `input_BAM`: Specify the full path to the final trimmed FASTQ file.
-  `LongReadSum_path` : Specify the full path to the LongReadSum software. 
3. **Expected Output:** See below.

| Output File Names      | Description |
| -------------- | ----------- |
| `bam_summary.txt`          | Summary text file containing statistics of the BAM file.  |
| `img`          | Directory containing images in the HTML output. |
| `st_bam_statistics_dynamic.html`    | HTML file with dynamic statistics of the BAM file.  |
| `st_bam_statistics.html`            | HTML file with static statistics of the BAM file. |





## 5. Gene Fusion Detection (Part 1/3) with LongGF `05_LongGF.sh`
1. **Goal:** Gene Fusion detection with LongGF
2. **Example Command:** `bash 05_LongGF.sh sample_name output_dir input_BAM_by_name ref_gtf min_support`
- `sample_name`: Enter the name for the sample.
- `output_dir`:  Specify the full path to 05_LongGF output. 
- `input_BAM`: Specify the full path to the BAM file sorted by read names.
- `ref_gtf`: Specify the full path to the reference GTF file containing gene annotations.
- `min_support`: Minimum support for gene fusion detection.
  
3. **Expected Output:** See below.

| Output File Names      |  | | | 
| -------------- | ----------- | ----------- | ----------- | 
| `longgf.100-25-100.txt`      | Min Overlap: 100bp | Bin Size: 25 | Min Map Length: 100bp |
| `longgf.100-10-100.txt`      | Min Overlap: 100bp | Bin Size: 10 | Min Map Length: 100bp |
| `longgf.75-25-75.txt`        | Min Overlap: 75bp | Bin Size: 25 | Min Map Length: 75bp  |
| `longgf.75-10-75.txt`        | Min Overlap: 75bp | Bin Size: 10 | Min Map Length: 75bp   |
| `longgf.50-25-50.txt`        | Min Overlap: 50bp | Bin Size: 25 | Min Map Length: 50bp   |
| `longgf.50-10-50.txt`        | Min Overlap: 50bp | Bin Size: 10 | Min Map Length: 50bp   |
| `longgf.25-25-25.txt`        | Min Overlap: 25bp | Bin Size: 25 | Min Map Length: 25bp   |
| `longgf.25-10-25.txt`        | Min Overlap: 25bp | Bin Size: 10 | Min Map Length: 25bp   |
| `longgf.10-25-10.txt`        | Min Overlap: 10bp | Bin Size: 25 | Min Map Length: 10bp   |
| `longgf.10-10-10.txt`        | Min Overlap: 10bp | Bin Size: 10 | Min Map Length: 10bp   |




## 5. Gene Fusion Detection (Part 2/3) with JAFFA `05_JAFFA.sh`
1. **Goal:** Gene Fusion detection with JAFFA
2. **Example Command:** `bash 05_JAFFA.sh sample_name input_fastq_gz output_dir JAFFA_bpipe JAFFAL_groovy_file`
- `sample_name`: Enter the name for the sample.
- `input_fastq_gz`: Specify the full path to the final trimmed and compressed FASTQ file.
- `output_dir`:  Specify the full path to 05_JAFFA output. 
- `JAFFA_bpipe`: Specify the path to the JAFFA bpipe executable (path/to/JAFFA/bpipe).
- `JAFFAL_groovy_file`: Specify the path to the JAFFA Long Read pipeline script (JAFFAL.groovy).
  
3. **Expected Output:** See below.

| Output File Names      | Description |
| -------------- | ----------- |
| `checks`                     | Directory containing checks by JAFFA.     |
| `commandlog.txt`             | Log file.  |
| `jaffa_results.csv`          | CSV file containing JAFFA fusion detection results. |
| `jaffa_results.fasta`        | FASTA file containing fusion sequences identified by JAFFA. |
| `sample_name_pass_reads_final_trimmed.fastq` | Directory containing various output files generated by JAFFA for the input FASTQ file. |






## 5. Gene Fusion Detection (Part 3/3) with FusionSeeker `05_FusionSeeker.sh`
1. **Goal:** Gene Fusion detection with FusionSeeker
2. **Example Command:** `bash 05_FusionSeeker.sh sample_name output_dir bam_input ref_gtf ref_genome`
- `sample_name`: Enter the name for the sample.
- - `output_dir`:  Specify the full path to 05_FusionSeeker output. 
- `bam_input`: Specify the full path to the input BAM file sorted by position.
- `ref_gtf`: Specify the full path to the reference GTF file containing gene annotations.
-  `ref_genome`: Specify the full path to the hg38 reference genome file (in FASTA format).
  
3. **Expected Output:** See below.

| Output File Names      | Description |
| -------------- | ----------- |
| `clustered_candidate.txt`                                | Clustered gene fusion candidates identified by FusionSeeker. |
| `confident_genefusion_transcript_sequence.fa`            | FASTA file containing sequences of confident gene fusion transcripts identified by FusionSeeker. |
| `confident_genefusion.txt`                               | File listing confident gene fusion events identified by FusionSeeker. |
| `log.txt`                                                | Log file. |
| `rawsignal.txt`                                          | File containing raw signal information from the FusionSeeker analysis. |




## 6. Gene Fusion Filtering Criteria `06_GF_criteria.py`
1. **Goal:** Combine the Gene Fusion results from the 3 programs together and filter based on:
- Genes of Interest
- Number of Supporting Reads
- Cross Reference to known Fusion Literature
  
2. **Script Customization Parameters:**
- 
  
3. **Expected Output:** See below.


