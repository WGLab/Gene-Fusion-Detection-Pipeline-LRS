#!/bin/bash

# Gene Fusion Detection 
# Purpose: Gene Fusion Detection with FusionSeeker 


sample_name=$1
output_dir=$2  #path/to/"05_FusionSeeker"


# Navigate to the output directory
cd "${output_dir}" || exit

mkdir -p "${sample_name}"
cd "${sample_name}" || exit


bam_input=$3   #<path/to/${sample_name}-guppy6.sorted_position.bam>

ref_gtf=$4  #<path/to/hg38/reference/gtf/annotation># gtf file

ref_genome=$5  # <path/to/hg38.fa>



#FusionSeeker
fusionseeker --bam $bam_input --datatype nanopore --gtf $ref_gtf $ref_genome



