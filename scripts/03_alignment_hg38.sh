#!/bin/bash

# Script to align the passed and Illumina adapter trimmed reads to the hg38 reference 
# Purpose: Splice aware alignment to hg38 reference using Minimap2 

# ALIGNMENT

sample_name=$1
input_fastq=$2  #<path/to/${sample_name}_pass_reads_final_trimmed.fastq>


output_dir=$3 #path/to/"03_alignment_hg38"


# reference files
REF_GENOME=$4 #<path/to/hg38.fa>
REF_SPLICE=$5 #<path/to/gencode.v44.annotation.bed> # bed file



# Navigate to the output directory

# Naviagate to the output directory and generate output directory 
cd "${output_dir}" || exit 
mkdir -p "${sample_name}" 

cd "${sample_name}" || exit


# Perform alignment using minimap2
minimap2 -ax splice --junc-bed ${REF_SPLICE} -uf -k14 -t 24 ${REF_GENOME} ${input_fastq} -o "${sample_name}-guppy6.sam"


# Convert SAM to BAM
samtools view -Sb -o "${sample_name}-guppy6.bam" "${sample_name}-guppy6.sam"

# Sort BAM file by position
samtools sort "${sample_name}-guppy6.bam" -o "${sample_name}-guppy6.sorted_position.bam"

# Index the sorted BAM file
samtools index "${sample_name}-guppy6.sorted_position.bam"

# Sort BAM by NAME for LongGF 
samtools sort -n "${sample_name}-guppy6.bam" -o "${sample_name}-guppy6.sorted_name.bam"

# Display completion message
echo "Alignment completed successfully for ${sample_name}"



