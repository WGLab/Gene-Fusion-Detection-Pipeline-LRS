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



## GF Detection Pipeline 

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
bash -o 01_redo_basecalling.out 01_redo_basecalling.sh example_sample14 /path/to/subset_example_sample14_all 01_redo_basecalling /path/to/guppy_basecaller config_file.cfg num_callers gpu_runners_per_device cpu_threads_per_caller device_cuda
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






