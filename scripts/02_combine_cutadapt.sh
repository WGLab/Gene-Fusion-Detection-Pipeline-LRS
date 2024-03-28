#!/bin/bash

# Script to first combine the passed reads into a single file as input to Cutadapt for Illumina Adapter Trimming  
# Purpose: Combine the passed reads after basecalling into a single file and Trim Illumina sequencing adapters from the reads 

sample_name=$1 
pass_input_dir=$2  #<path/to/re-basecalled/pass/folder>

output_dir=$3  #path/to/"02_combine_cutadapt"


# Cutadapt adapters 
adapter3=$4
adapter5=$5

rev_adapter3=$6 
rev_adapter5=$7

error_rate_trim1=$8 
error_rate_trim2=$9

num_occurances_trim1=${10}
num_occurances_trim2=${11}



# Naviagate to the output directory and generate output directory 
cd "${output_dir}" || exit 
mkdir -p "${sample_name}" 

cd "${sample_name}" || exit


# Combine files
cat $pass_input_dir/f*.fastq > ${sample_name}_pass_reads.fastq

CUTADAPT_INPUT="${output_dir}/${sample_name}/${sample_name}_pass_reads.fastq"


# Cutadapt Adapter Trimming 
# PART 1 - trims adapters around region of interest - and outputs a trimmed and untrimmed file
# PART 2 - trims all adapters around that region of interest (will trim everything left and right of the adapter)
# doing 2 separate trims to account for the adapters stuck in the middle after the first initial trim at the edges

#NOTE: FASTQ FILE WITH CUTADAPT TRIMMING IS WITHIN THE SAME FOLDER AS THE CONCATENATED FASTQ FILE


##trim adapters if they 5' adapters and 3' adapters flank the sequence of interest
##output trimmed and untrimmed fastq
cutadapt -g "$adapter5...$adapter3" -g "$rev_adapter5...$rev_adapter3" \
        -e "$error_rate_trim1" -n "$num_occurances_trim1" \
        --untrimmed-output ${CUTADAPT_INPUT/.fastq/_untrimmed1.fastq} \
        -o ${CUTADAPT_INPUT/.fastq/_trimmed1.fastq} $CUTADAPT_INPUT >> "cutadapt_report_linked_adapters.txt"


##trim adapters regardless of if they flank the region of interest
##uses the untrimmed output of the first round of cutadapt as input
cutadapt -g $adapter5 -g $rev_adapter5 -a $adapter3 -a $rev_adapter3 \
        -e "$error_rate_trim2" -n "$num_occurances_trim2" \
        --untrimmed-output ${CUTADAPT_INPUT/.fastq/_untrimmed.fastq} \
        -o ${CUTADAPT_INPUT/.fastq/_trimmed2.fastq} ${CUTADAPT_INPUT/.fastq/_untrimmed1.fastq} >> "cutadapt_report_single_adapters.txt"

cat ${CUTADAPT_INPUT/.fastq/_trimmed1.fastq} ${CUTADAPT_INPUT/.fastq/_trimmed2.fastq} >> ${CUTADAPT_INPUT/.fastq/_final_trimmed.fastq}

# gzip final trimmed file
gzip -k ${CUTADAPT_INPUT/.fastq/_final_trimmed.fastq}

