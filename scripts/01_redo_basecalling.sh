#!/bin/bash

# Script to Redo Basecalling from long-reads sequenced and basecalled with ONT's integrated MinKNOW basecaller
# Purpose: Re-running basecalling to help increase the accuracy of the sequenced reads 

sample_name=$1
input_dir=$2
output_dir=$3 # path/to/"01_redo_basecalling"


basecaller_input=$4  #path/to/guppy6/basecaller_model 
config_file=$5 
num_callers=$6
gpu_runners_per_device=$7 
cpu_threads_per_caller=$8
device_cuda=$9




# Naviagate to the output directory and generate output directory
cd "${output_dir}" || exit
mkdir -p "${sample_name}"

cd "${sample_name}" || exit



# Redo Basecalling 

${basecaller_input} -i ${input_dir} \
-s "${output_dir}/${sample_name}" \
-c ${config_file}> \
--recursive --num_callers ${num_callers} --gpu_runners_per_device ${gpu_runners_per_device} --cpu_threads_per_caller ${cpu_threads_per_caller} --device cuda:${device_cuda} 


echo "Redo basecalling is completed."

