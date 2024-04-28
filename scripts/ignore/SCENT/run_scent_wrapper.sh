#!/bin/bash

dir=/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/info

projName=tmparm3
celltype_to_use=HSCs_all
num=${celltype_to_use}_nums
downsample=FALSE
use_interaction=TRUE
line=185

while read -r line; do
sbatch --account=smontgom --partition=batch --time=1-1:00:00 --mem=4G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21_multiome/scripts/andrew/SCENT/run_scent.sh $projName $celltype_to_use $line $downsample $use_interaction
echo $line
done < $dir/$num

# for line in 194 201 202; do
# sbatch --account=smontgom --partition=batch --time=1-1:00:00 --mem=4G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21_multiome/scripts/andrew/SCENT/run_scent.sh $projName $celltype_to_use $line $downsample $use_interaction
# done