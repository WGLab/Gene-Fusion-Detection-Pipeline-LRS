# Example using a subset of Sample #14 

## Background
The example dataset is a subset of reads derived from Sample #14 and consists of a total of 3,064 reads previously aligned to gene breakpoints of the expected gene fusion (GF), _KIAA1549_::_BRAF_. This included reads that aligned to the reported breakpoints of this expected GF, which includes reads that support the expected fusion and those that do not. This example offers a detailed examination and explanation of the pipeline usage including both reads that do and do not necessarily support the expected GF. 

## Input Processing 
These reads were first extracted from the previously aligned bam files in the original analysis, using `samtools view` based on the region surrounding the breakpoint by 100 bps on either side to be used in the example, then subsetted from the original fast5 files using `fast5_subset` and then compressed into a single file named, `subset_example_sample14_all.tar.gz`. 

Command to make a new folder to store the contents of after decompressing the file of input reads: 
```
mkdir subset_example_sample14_all && tar -xzvf subset_example_sample14_all.tar.gz -C subset_example_sample14_all
```

When you run this command, the output will be a list of files in the newly created `subset_example_sample14_all` folder, as seen below.
```
batch0.fast5   batch12.fast5  batch15.fast5  batch18.fast5  batch20.fast5  batch23.fast5  batch26.fast5  batch29.fast5  batch4.fast5  batch7.fast5
batch10.fast5  batch13.fast5  batch16.fast5  batch19.fast5  batch21.fast5  batch24.fast5  batch27.fast5  batch2.fast5   batch5.fast5  batch8.fast5
batch11.fast5  batch14.fast5  batch17.fast5  batch1.fast5   batch22.fast5  batch25.fast5  batch28.fast5  batch3.fast5   batch6.fast5  batch9.fast5
```
Note the path to this folder for later.



# GF Detection Pipeline 

## 0. Initialization of Directories `00_dir_setup.sh`
In the folder where the downloaded scripts are stored, run the `00_dir_setup.sh` script to perform the initialization of the necessary directories, as seen below. 

**Example Command:** `bash 00_dir_setup.sh`

**Expected Output:**
```
01_redo_basecalling  02_combine_cutadapt  03_alignment_hg38  04_longreadsum  05_FusionSeeker  05_JAFFA  05_LongGF  06_GF_criteria
```


## 1. Redo Basecalling with Guppy6 `01_redo_basecalling.sh`

**Base Command:** `bash 01_redo_basecalling.sh sample_name input_dir output_dir basecaller_input config_file num_callers gpu_runners_per_device cpu_threads_per_caller device_cuda`
- `sample_name`: Enter the name for the sample.
- `input_dir`: Specify the path to the input directory.
- `output_dir`: Specify the full path to 01_redo_basecalling output. 
- `basecaller_input`: Specify the path to the basecaller model. 
- `config_file`: Specify the path to the configuration file name.
- Additional HPC Machine Configuration (Adjust based on your HPC machine specifications):
   - `num_callers`, `gpu_runners_per_device`, `cpu_threads_per_caller`, and `device_cuda`
 
**Example Command:** 
```
bash 01_redo_basecalling.sh example_sample14 /path/to/subset_example_sample14_all 01_redo_basecalling /path/to/guppy_basecaller config_file.cfg num_callers gpu_runners_per_device cpu_threads_per_caller device_cuda
```

**Expected Output:** Within the '01_redo_basecalling' folder 
```
guppy_basecaller_log-2024-03-28_13-09-44.log  pass  sequencing_summary.txt  sequencing_telemetry.js
```
There is one fastq file within the pass folder named, `fastq_runid_d1615b727f9b8bba03d838cc3386028b4a8b9ab0_0_0.fastq`, which consists of the rebasecalled reads in a single fastq file. 






## 2. Combine and Cutadapt `02_combine_cutadapt.sh`

**Base Command:** `bash 02_combine_cutadapt.sh sample_name pass_input_dir output_dir adapter3 adapter5 rev_adapter3 rev_adapter5 error_rate_trim1 error_rate_trim2 num_occurances_trim1 num_occurances_trim2`
- `sample_name`: Enter the desired name for the sample.
- `pass_input_dir`: Specify the full path to the "pass" folder of the re-basecalled samples.
- `output_dir`: Specify the full path to 02_combine_cutadapt output. 
- `adapter3`, `adapter5`, `rev_adapter3`, and `rev_adapter5`: Forward and Reverse adapter sequences for trimming.
- `error_rate_trim1` and `error_rate_trim2`: Cutadapt error rate denoted for both trimming steps. 
- `num_occurances_trim1` and `num_occurances_trim2`: Reads with fewer than the specified occurrences of a specified adapter were trimmed, denoted for both trimming steps. 

**Example Command:** 
```
bash 02_combine_cutadapt.sh example_sample14 /path/to/01_redo_basecalling/example_sample14/pass /path/to/02_combine_cutadapt AGATCGGAAGAGCACACGTCTGAACTCCAGTCA ACACTCTTTCCCTACACGACGCTCTTCCGATCT AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT 0.3 0.3 10 20
```

**Expected Output:** Within the '02_combine_cutadapt' folder 
```
cutadapt_report_linked_adapters.txt  example_sample14_pass_reads_final_trimmed.fastq     example_sample14_pass_reads_trimmed2.fastq
cutadapt_report_single_adapters.txt  example_sample14_pass_reads_final_trimmed.fastq.gz  example_sample14_pass_reads_untrimmed1.fastq
example_sample14_pass_reads.fastq    example_sample14_pass_reads_trimmed1.fastq          example_sample14_pass_reads_untrimmed.fastq
```


### Example expected output from: `cutadapt_report_linked_adapters.txt`
```
=== Summary ===

Total reads processed:                   2,948
Reads with adapters:                     2,133 (72.4%)

== Read fate breakdown ==
Reads discarded as untrimmed:              815 (27.6%)
Reads written (passing filters):         2,133 (72.4%)

Total basepairs processed:     1,147,481 bp
Total written (filtered):        541,240 bp (47.2%)
```

### Example expected output from: `cutadapt_report_single_adapters.txt`
```
=== Summary ===

Total reads processed:                     815
Reads with adapters:                       815 (100.0%)

== Read fate breakdown ==
Reads discarded as untrimmed:                0 (0.0%)
Reads written (passing filters):           815 (100.0%)

Total basepairs processed:       304,568 bp
Total written (filtered):        220,695 bp (72.5%)
```



## 3. Alignment to the hg38 Reference `03_alignment_hg38.sh`

**Base Command:** `bash 03_alignment_hg38.sh sample_name input_fastq output_dir REF_GENOME REF_SPLICE`
- `sample_name`: Enter the name for the sample.
- `input_fastq`: Specify the full path to the final trimmed FASTQ file.
- `output_dir`: Specify the full path to 03_alignment_hg38 output. 
- `REF_GENOME`: Specify the full path to the hg38 reference genome file (in FASTA format).
- `REF_SPLICE`: Specify the full path to the BED file containing splice site annotations for hg38.

**Example Command:** 
```
bash 03_alignment_hg38.sh example_sample14 /path/to/02_combine_cutadapt/example_sample14/example_sample14_pass_reads_final_trimmed.fastq  /path/to/03_alignment_hg38 /path/to/ref_genome_hg38/hg38.fa /path/to/gencode.v44.annotation.bed
```


**Expected Output:** Within the '03_alignment_hg38' folder 
```
example_sample14-guppy6.bam  example_sample14-guppy6.sorted_name.bam      example_sample14-guppy6.sorted_position.bam.bai
example_sample14-guppy6.sam  example_sample14-guppy6.sorted_position.bam
```




## 4. Quality Check Analysis with LongReadSum  `04_longreadsum.sh`
**Example Command:** `bash 04_longreadsum.sh sample_name output_dir input_BAM long_read_sum_path`
- `sample_name`: Enter the name for the sample.
- `output_dir`:  Specify the full path to 04_longreadsum output. 
- `input_BAM`: Specify the full path to the aligned BAM file.
-  `LongReadSum_path` : Specify the full path to the LongReadSum software. 

**Example Command:** 
```
bash 04_longreadsum.sh example_sample14 /path/to/04_longreadsum /path/to/03_alignment_hg38/example_sample14/example_sample14-guppy6.bam /path/to/LongReadSum
```

**Expected Output:** Within the '04_longreadsum' folder 
```
bam_summary.txt  img  st_bam_statistics_dynamic.html  st_bam_statistics.html
```

### Example expected output from: `bam_summary.txt`
```
Longest read length     1527
N50 read length 280
Mean read length        258.46
Median read length      229
GC%     47.71

Total number of mapped reads    2948
Total number of mapped bases    761935
Longest mapped read length      1527
N50 mapped read length  280
Mean mapped read length 258.46
Median mapped read length       229
GC%     47.71

Total number of primary alignments      2948
Total number of secondary alignments    65
Total number of supplementary alignments        176
Total number of reads with secondary alignments 59
Total number of reads with supplementary alignments     166
Total number of reads with both secondary and supplementary alignments  16
Total number of reads with forward alignments   1671
Total number of reads with reverse alignments   1277
Total number of reverse alignment       1277

Total number of matched bases   667869
Total number of mismatched bases        39937
Total number of insertions      13855
Total number of deletions       13039
Total number of soft clipped bases      168149
```




## 5. Gene Fusion Detection (Part 1/3) with LongGF `05_LongGF.sh`
**Example Command:** `bash 05_LongGF.sh sample_name output_dir input_BAM_by_name ref_gtf min_support`
- `sample_name`: Enter the name for the sample.
- `output_dir`:  Specify the full path to 05_LongGF output. 
- `input_BAM`: Specify the full path to the BAM file sorted by read names.
- `ref_gtf`: Specify the full path to the reference GTF file containing gene annotations.
- `min_support`: Minimum support for gene fusion detection.


**Example Command:** 
```
bash 05_LongGF.sh example_sample14 /path/to/05_LongGF /path/to/03_alignment_hg38/example_sample14/example_sample14-guppy6.sorted_name.bam /path/to/gencode.v44.annotation.gtf 1
```

**Expected Output:** Within the '04_longreadsum' folder 
```
longgf.100-10-100.txt  longgf.10-10-10.txt  longgf.25-10-25.txt  longgf.50-10-50.txt  longgf.75-10-75.txt
longgf.100-25-100.txt  longgf.10-25-10.txt  longgf.25-25-25.txt  longgf.50-25-50.txt  longgf.75-25-75.txt
```


**Example Contents of `longgf.100-10-100.txt`:**
```
GF      KIAA1549:BRAF 10 10 supporting reads=10/10 138867973:9;138867977:1 140787583:10 chr7:138867973 10/0:24:34 chr7:140787583 0/10:24:185
        138867973(-chr7:138867973-138869738/03dec4eb-8641-40bc-9c86-61c02705affa:18-370)1 140787583(-chr7:140781667-140787583/370-568)1
        138867973(+chr7:138867973-138869586/2e24c422-39dd-4b6f-9a9b-ccb98b3cb3c4:201-400)1 140787583(+chr7:140781663-140787583/0-201)1
        138867973(-chr7:138867973-138869575/43cd9a89-9470-47f0-a43d-bdf81691fb51:21-215)1 140787583(-chr7:140783102-140787583/215-419)1
        138867973(+chr7:138867973-138868079/47d3a423-8a28-47bd-a50e-381ec6be022f:199-305)1 140787583(+chr7:140781663-140787583/0-199)1
        138867973(+chr7:138867973-138869693/7e56e365-97f1-46fa-bb1c-3ccb85045ea8:198-514)1 140787583(+chr7:140781663-140787583/0-198)1
        138867977(-chr7:138867977-138868129/b0c890e6-de00-488c-b208-f943cefa23fe:32-185)1 140787583(-chr7:140783084-140787583/188-299)1
        138867973(-chr7:138867973-138871335/b0e9acd9-8e33-45f9-91d9-cefcef3063f3:17-565)1 140787583(-chr7:140781663-140787583/565-772)1
        138867973(-chr7:138867973-138868081/c01d20d4-b3fc-47a3-b6ca-1f33fa766df4:21-128)1 140787583(-chr7:140781663-140787583/128-329)1
        138867973(+chr7:138867973-138871212/c63aafef-89ab-4c0f-aa5c-a0025b8e5957:198-624)1 140787583(+chr7:140781663-140787583/0-198)1
        138867973(+chr7:138867973-138868129/ffd7691d-da58-4bc5-a591-36a250797cae:201-358)1 140787583(+chr7:140781663-140787583/0-201)1
SumGF   KIAA1549:BRAF 10 chr7:138867973 chr7:140787583
GF      BRAF:BCR 1 1 supporting reads=1/1 140783032:1 23310417:1 chr7:140783032 1/0:2:185 chr22:23310417 0/1:2:1
        140783032(+chr7:140783032-140794348/6272f2c5-e8d8-4902-9ddd-74cbc256eca5:227-429)1 23310417(+chr22:23295078-23310417/18-243)1
SumGF   BRAF:BCR 1 chr7:140783032 chr22:23310417
```







## 5. Gene Fusion Detection (Part 2/3) with JAFFA `05_JAFFA.sh`
**Example Command:** `bash 05_JAFFA.sh sample_name input_fastq_gz output_dir JAFFA_bpipe JAFFAL_groovy_file`
- `sample_name`: Enter the name for the sample.
- `input_fastq_gz`: Specify the full path to the final trimmed and compressed FASTQ file.
- `output_dir`:  Specify the full path to 05_JAFFA output. 
- `JAFFA_bpipe`: Specify the path to the JAFFA bpipe executable (path/to/JAFFA/bpipe).
- `JAFFAL_groovy_file`: Specify the path to the JAFFA Long Read pipeline script (JAFFAL.groovy).



**Example Command:** 
```
bash 05_JAFFA.sh example_sample14 /path/to/02_combine_cutadapt/example_sample14/example_sample14_pass_reads_final_trimmed.fastq.gz /path/to/05_JAFFA /path/to/JAFFA-version-2.3/tools/bpipe-0.9.9.2/bin/bpipe /path/to/JAFFA-version-2.3/JAFFAL.groovy
```

**Expected Output:** Within the '04_longreadsum' folder 
```
checks  commandlog.txt  example_sample14_pass_reads_final_trimmed.fastq  jaffa_results.csv  jaffa_results.fasta
```

**Example Contents of `jaffa_results.csv`:**
```
sample,fusion genes,chrom1,base1,strand1,chrom2,base2,strand2,gap (kb),spanning pairs,spanning reads,inframe,aligns,rearrangement,contig,contig break,classification,known
example_sample14_pass_reads_final_trimmed.fastq,KIAA1549:BRAF,chr7,138867975,-,chr7,140787584,-,1919.608,0,10,TRUE,TRUE,TRUE,03dec4eb-8641-40bc-9c86-61c02705affa,370,HighConfidence,Yes
example_sample14_pass_reads_final_trimmed.fastq,KIAA1549:BRAF,chr7,138867975,-,chr7,140783192,-,1915.217,0,33,NA,FALSE,TRUE,0ee394f5-e4a1-4c43-8293-edc332d6b48b,90,LowConfidence,Yes
```





## 5. Gene Fusion Detection (Part 3/3) with FusionSeeker `05_FusionSeeker.sh`
**Example Command:** `bash 05_FusionSeeker.sh sample_name output_dir bam_input ref_gtf ref_genome`
- `sample_name`: Enter the name for the sample.
- `output_dir`:  Specify the full path to 05_FusionSeeker output. 
- `bam_input`: Specify the full path to the input BAM file sorted by position.
- `ref_gtf`: Specify the full path to the reference GTF file containing gene annotations.
-  `ref_genome`: Specify the full path to the hg38 reference genome file (in FASTA format).

**Example Command:** 
```
bash 05_FusionSeeker.sh example_sample14 /path/to/05_FusionSeeker /path/to/03_alignment_hg38/example_sample14/example_sample14-guppy6.sorted_position.bam /path/to/gencode.v44.annotation.gtf /path/to/hg38.fa
```

**Expected Output:** Within the '04_longreadsum' folder 
```
clustered_candidate.txt  confident_genefusion.txt  log.txt  poa_workspace  raw_signal  rawsignal.txt
```

**Example Contents of `confident_genefusion.txt`:**
```
KIAA1549        BRAF    11      chr7    138867974       chr7    140787582       GF01    ffd7691d-da58-4bc5-a591-36a250797cae,03dec4eb-8641-40bc-9c86-61c02705affa,1de29aef-bd98-47a1-b866-49843ddf3fdb,2e24c422-39dd-4b6f-9a9b-ccb98b3cb3c4,43cd9a89-9470-47f0-a43d-bdf81691fb51,47d3a423-8a28-47bd-a50e-381ec6be022f,7e56e365-97f1-46fa-bb1c-3ccb85045ea8,b0c890e6-de00-488c-b208-f943cefa23fe,b0e9acd9-8e33-45f9-91d9-cefcef3063f3,c01d20d4-b3fc-47a3-b6ca-1f33fa766df4,c63aafef-89ab-4c0f-aa5c-a0025b8e5957
```







## 6. Gene Fusion Filtering Criteria `06_GF_criteria.py`
**Example Command:** `python 06_GF_criteria.py ref_gtf_annotation sample_name artifact_filter_path input_longgf input_jaffa input_fusion_seek bam_path genes_of_interest num_supporting_reads mitelman_db_gfs chimer_db_gfs cosmic_db_gfs > 06_criteria.out`
- `ref_gtf_annotation`: path/to/gencode.v44.annotation.bed
- `sample_name`: Enter the name of the sample.
- `artifact_filter_path`: Specify full path to 06_GF_criteria output folder.
- `input_longgf`: Specify the full path to the sample_name 05_LongGF folder containing the different outputs of the LongGF parameters.
- `input_jaffa`: Specify the full path to 05_JAFFA folder. 
- `input_fusion_seek`: Specify the full path to 05_FusionSeeker folder. 
- `bam_path`: Specify the full path to the sorted by position BAM file.
- `genes_of_interest`: Specify the full path to the genes of interest file to prioritize the results.
- `num_supporting_reads`: Specify the minimum number of supporting reads threshold.
- `mitelman_db_gfs`: Specify the path to mitelman_db_fusion_names.txt
- `chimer_db_gfs`: Specify the path to the ChimerDB4.0_fusions.txt
- `cosmic_db_gfs`: Specify the full path to the cosmic_fusion_genes.tsv


**Example Command:** 
```
python 06_GF_criteria.py /path/to/gencode.v44.annotation.bed example_sample14 /path/to/06_GF_criteria/ /path/to/05_LongGF/ /path/to/05_JAFFA/ /path/to/05_FusionSeeker/ /path/to/03_alignment_hg38/example_sample14/example_sample14-guppy6.sorted_position.bam /path/to/CHOP_Fusion_Panel_genes.txt 2 /path/to/fusion_names.txt /path/to/ChimerDB4.0_fusions.txt /path/to/cosmic_fusion_genes.tsv > 06_criteria.out 
```
- note the backslashes are necessary in the input and output paths specified in the command


**Expected Output:** Within the '04_longreadsum' folder 
```
all_fusions_genes_of_interest_filtered.tsv  all_fusions_num_support_reads_filtered.tsv  example_sample14_all_longgf.tsv    example_sample14_unique_longgf.tsv  literature_cross_referenced_output.tsv
all_fusions_genes_of_interest.tsv           example_sample14_all_FusionSeeker.tsv       example_sample14_fusionseeker.tsv  final_output_no_readIDs.tsv         longgf_unique_orientation_check.tsv
all_fusions_LGF_J_FS.tsv                    example_sample14_all_JAFFA.tsv              example_sample14_jaffa.tsv         final_output.tsv                    longgf_unique_strandness.tsv
```

**Example Contents of `final_output.tsv` Folder:**
```
Fusion Name     Gene 1 Symbol   Gene 1 Breakpoint       Gene 2 Symbol   Gene 2 Breakpoint       Supporting Reads        Read IDs        Literature Cross Reference      GF Program
KIAA1549:BRAF   KIAA1549        chr7:138867973  BRAF    chr7:140787583  40      03dec4eb-8641-40bc-9c86-61c02705affa,08c6f667-8e89-4e10-b013-a0a0c0f43744,0ee394f5-e4a1-4c43-8293-edc332d6b48b,1b37b834-49f5-4055-9e35-5c495a241751,1f88d9b5-4c28-48d4-b8ea-c80b3ea2f946,249369e5-ea14-4bb7-acbe-c2d8e11ae2e1,28c21dcd-1b00-4207-b5eb-f6b43e5cedbd,2e24c422-39dd-4b6f-9a9b-ccb98b3cb3c4,2fb319a5-bd11-4c04-b480-8f120bdb28b9,3ba42710-a6a2-4fcd-9c92-d44b148fed6a,3fc80036-1bfc-49ae-bd63-d424e6456495,43cd9a89-9470-47f0-a43d-bdf81691fb51,47d3a423-8a28-47bd-a50e-381ec6be022f,4f54aea0-6bbc-45ae-b0ad-6fa976ecfc1e,5156b766-fe9e-490c-aee6-697d694e85c5,5a719cda-b909-4212-8033-26d264d14389,618a60ec-5650-4c24-914e-28228dd6315a,63eb7086-a2f1-482c-b048-6856ecd2b0ed,6b6b8f23-f234-4adc-98cb-a355736feb72,7b950612-a91e-4afb-9eb3-73659d0285f4,7e56e365-97f1-46fa-bb1c-3ccb85045ea8,8e91e67e-98a7-4408-877a-cc6c7620ebeb,9b2a9942-2c38-432a-9263-36d89dde415b,9c4db01f-2f85-41f0-a843-4babf645078d,ae44809d-9074-4865-96ab-c8db51b3015c,af3a6867-8136-4a21-b371-29704d3d3128,b0c890e6-de00-488c-b208-f943cefa23fe,b0e9acd9-8e33-45f9-91d9-cefcef3063f3,b4ce4c4c-c70a-4714-a949-3d25a99fb58e,bc693a9b-48ee-44d0-be34-8a6825246188,c01d20d4-b3fc-47a3-b6ca-1f33fa766df4,c14fb703-6d11-4151-9cdb-6eda70c4a85e,c63aafef-89ab-4c0f-aa5c-a0025b8e5957,c732d8da-cb72-4cb6-81ec-09fb136dbce9,cb2b8543-748b-41f4-be41-96ca3fffca68,cea878af-060f-40b2-8bec-ff2812dcee98,eee6ee9d-0d87-49ab-ae55-fbfbd03b6fbc,f64b8865-79f2-4478-9128-2f88075b4622,f88b6c33-4915-4d56-a6ca-d07d24453c21,ffd7691d-da58-4bc5-a591-36a250797cae  KIAA1549:BRAF, KIAA1549:BRAF    LongGF
KIAA1549:BRAF   KIAA1549        chr7:138867975  BRAF    chr7:140783192  33      0ee394f5-e4a1-4c43-8293-edc332d6b48b    KIAA1549:BRAF, KIAA1549:BRAF    JAFFA
KIAA1549:BRAF   KIAA1549        chr7:138867974  BRAF    chr7:140787582  11      GF01    KIAA1549:BRAF, KIAA1549:BRAF    FusionSeeker
KIAA1549:BRAF   KIAA1549        chr7:138867975  BRAF    chr7:140787584  10      03dec4eb-8641-40bc-9c86-61c02705affa    KIAA1549:BRAF, KIAA1549:BRAF    JAFFA
```


### From review of the final_output.tsv, we are able to see that the expected causal GF in this example dataset is **_KIAA1549_::_BRAF_**. 
