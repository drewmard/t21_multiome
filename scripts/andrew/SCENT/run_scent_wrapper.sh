#!/bin/bash

dir=/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/info

projName=tmparm2
celltype_to_use=HSCs_H
num=${celltype_to_use}_nums

while read -r line; do
sbatch --account=smontgom --partition=batch --time=1-1:00:00 --mem=4G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/t21_multiome/scripts/andrew/SCENT/run_scent.sh $projName $celltype_to_use $line
echo $line
done < $dir/$num