# srun --account=default --partition=interactive --time=24:00:00 --mem=128G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash

library(Seurat)
library(Signac)
f="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds"
df = readRDS(f)

SplitFragments(subset(df,disease=="H" & subclust_v6=="HSCs"),group.by="subclust_v6",outdir = "/oak/stanford/groups/smontgom/amarder/tmp",assay = "ATAC",file.suffix='disomy')
SplitFragments(subset(df,disease=="T21" & subclust_v6=="HSCs"),group.by="subclust_v6",outdir = "/oak/stanford/groups/smontgom/amarder/tmp",assay = "ATAC",file.suffix='ts21')
SplitFragments(subset(df,subclust_v6=="HSCs"),group.by="subclust_v6",outdir = "/oak/stanford/groups/smontgom/amarder/tmp",assay = "ATAC",file.suffix='all')

# # run the next in bash:
# # 
# awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+";  print $1,mid,$3,"N",1000,"-"}' /oak/stanford/groups/smontgom/amarder/tmp/HSCsall.bed | sort -k 1,1V -k2,2n | bgzip > /oak/stanford/groups/smontgom/amarder/tmp/HSCsall.tagAlign.gz
# 
# awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+";  print $1,mid,$3,"N",1000,"-"}' /oak/stanford/groups/smontgom/amarder/tmp/HSCsts21.bed | sort -k 1,1V -k2,2n | bgzip > /oak/stanford/groups/smontgom/amarder/tmp/HSCsts21.tagAlign.gz
#
# # awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+";  print $1,mid,$3,"N",1000,"-"}' /oak/stanford/groups/smontgom/amarder/tmp/HSCsdisomy.bed | sort -k 1,1V -k2,2n | bgzip > /oak/stanford/groups/smontgom/amarder/tmp/HSCsdisomy.tagAlign.gz
