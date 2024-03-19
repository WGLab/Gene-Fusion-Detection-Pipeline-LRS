#!/bin/bash

# Quality check analysis with LongReadSum 
# Purpose: Perform a quality check analysis of the mapping to hg38 reference with LongReadSum 

sample_name=$1 
output_dir=$2 #path/to/"04_longreadsum"

# Navigate to the output directory
cd "${output_dir}" || exit

mkdir -p "${sample_name}"
cd "${sample_name}" || exit


input_BAM=$3  #<path/to/${sample_name}-guppy6.sorted_position.bam>

long_read_sum_path=$4  # <path/to/software/LongReadSum>

#LongReadSum 
python ${long_read_sum_path} bam -i $input_BAM


