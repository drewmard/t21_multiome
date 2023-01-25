#!/bin/bash

# # sbatch --account=smontgom --partition=batch --time=1-1:00:00 --mem=72G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/RV/run_create_peak_RV_dataset.sh

module load R/4.1.2

projName=$1
celltype_to_use=$2
num=$3
downsample=$4

path_to_script=/oak/stanford/groups/smontgom/amarder/t21_multiome/scripts/andrew/SCENT/run_scent.R 
Rscript $path_to_script $projName $celltype_to_use $num $downsample

