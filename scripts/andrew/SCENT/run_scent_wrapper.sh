#!/bin/bash

dir=/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/info

projName=tmparm3
celltype_to_use=HSCs_all
num=${celltype_to_use}_nums
downsample=FALSE

while read -r line; do
sbatch --account=smontgom --partition=batch --time=1-1:00:00 --mem=4G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21_multiome/scripts/andrew/SCENT/run_scent.sh $projName $celltype_to_use $line $downsample
echo $line
done < $dir/$num
