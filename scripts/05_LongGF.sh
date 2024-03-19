#!/bin/bash

# Gene Fusion Detection 
# Purpose: Gene Fusion Detection with LongGF across 10 different parameters 

sample_name=$1

output_dir=$2


input_BAM_by_name=$3  #<path/to/${sample_name}-guppy6.sorted_name.bam>

ref_gtf=$4  #<path/to/hg38/reference/gtf/annotation># gtf file

min_support=$5



# Navigate to the output directory
cd "${output_dir}" || exit

mkdir -p "${sample_name}"
cd "${sample_name}" || exit



# LongGF 
# 100bp overlap
LongGF $input_BAM_by_name $ref_gtf 100 25 100 0 0 ${min_support} 0 > longgf.100-25-100.txt
LongGF $input_BAM_by_name $ref_gtf 100 10 100 0 0 ${min_support} 0 > longgf.100-10-100.txt

# 75bp overlap 
LongGF $input_BAM_by_name $ref_gtf 75 25 75 0 0 ${min_support} 0 > longgf.75-25-75.txt
LongGF $input_BAM_by_name $ref_gtf 75 10 75 0 0 ${min_support} 0 > longgf.75-10-75.txt

# 50bp overlap
LongGF $input_BAM_by_name $ref_gtf 50 25 50 0 0 ${min_support} 0 > longgf.50-25-50.txt
LongGF $input_BAM_by_name $ref_gtf 50 10 50 0 0 ${min_support} 0 > longgf.50-10-50.txt

# 25bp overlap
LongGF $input_BAM_by_name $ref_gtf 25 25 25 0 0 ${min_support} 0 > longgf.25-25-25.txt
LongGF $input_BAM_by_name $ref_gtf 25 10 25 0 0 ${min_support} 0 > longgf.25-10-25.txt

# 10bp overlap
LongGF $input_BAM_by_name $ref_gtf 10 25 10 0 0 ${min_support} 0 > longgf.10-25-10.txt
LongGF $input_BAM_by_name $ref_gtf 10 10 10 0 0 ${min_support} 0 > longgf.10-10-10.txt

# Display completion message
echo "LongGF detection completed successfully for ${sample_name}"


