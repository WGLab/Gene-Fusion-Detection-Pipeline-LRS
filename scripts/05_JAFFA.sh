#!/bin/bash

# Gene Fusion Detection 
# Purpose: Gene Fusion Detection with JAFFA 


sample_name=$1 
input_fastq_gz=$2  #<path/to/${sample_name}_pass_reads_final_trimmed.fastq.gz>

output_dir=$3  #path/to/"05_JAFFA"


JAFFA_bpipe=$4  #path/to/JAFFA/bpipe
JAFFAL_groovy_file=$5  #path/to/JAFFAL.groovy

# Navigate to the output directory
cd "${output_dir}" || exit

mkdir -p "${sample_name}"
cd "${sample_name}" || exit


# JAFFA long parameters
${JAFFA_bpipe} run ${JAFFAL_groovy_file} $input_fastq_gz


