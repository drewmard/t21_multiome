# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# module load R/4.1.2

# Code adapted from SCENT: https://github.com/immunogenomics/SCENT/blob/main/SCENT.R (original author: Saori Sakaue)

# library(biomaRt)
library(Seurat)
library(Signac)
library(data.table)
library(lme4)
library(stringr)
library(boot)
library(MASS)
library(Matrix)
library(parallel)
options(stringsAsFactors = F)

downsample=FALSE
source("/oak/stanford/groups/smontgom/amarder/t21_multiome/scripts/andrew/SCENT/scent_v4.R")

###############################################################################
###############################################################################
###############################################################################


projName="tmparm"
celltype_to_use = "HSCs_T21"
num=77

args = commandArgs(trailingOnly=TRUE)
projName = args[1]
celltype_to_use = args[2]
num = as.numeric(args[3])
if (length(args) > 3) {downsample=args[4]} else {downsample="FALSE"}
print(paste("run_scent.R",projName,celltype_to_use,num,downsample))
fDir = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input"
tmpDir = paste0("/oak/stanford/groups/smontgom/amarder/tmp/",projName)
f.atac_out = paste0(tmpDir,"/atac/",celltype_to_use,".atac.",num,".rds")
f.rna_out = paste0(tmpDir,"/rna/",celltype_to_use,".rna.",num,".rds")
f.meta_out = paste0(fDir,"/meta/meta.",celltype_to_use,".rds")
f.chunkinfo = paste0(tmpDir,"/chunkinfo_split/",celltype_to_use,".chunkinfo.",num,".txt")
output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out_split/",celltype_to_use,".",num,".txt")

atac.all = readRDS(f.atac_out)
mrna = readRDS(f.rna_out)
meta = readRDS(f.meta_out)
chunkinfo = fread(f.chunkinfo,data.table = F,stringsAsFactors = F)

if (downsample=="TRUE") {
  print("~ ~ ~ DOWNSAMPLING ~ ~ ~")
  set.seed(03191995)
  ind = sort(sample(1:3784,2431,replace = F))
  atac.all = atac.all[,ind]
  mrna = mrna[,ind]
  meta = meta[ind,]
  output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out_split/",celltype_to_use,".down.",num,".txt")
}

print(paste0("Running ",nrow(chunkinfo)," peak-gene connections..."))
system.time(out.sub <- lapply(1:nrow(chunkinfo),
                              create_input_and_run_SCENT,
                              run_bs=TRUE,
                              bootstrap_sig=TRUE))
# numThreads = detectCores() #/2
# system.time(out.sub <- mclapply(1:3,
#                                 create_input_and_run_SCENT,
#                                 run_bs=TRUE,
#                                 bootstrap_sig=TRUE,
#                                 mc.cores = min(nrow(chunkinfo),numThreads)
#                                 ))
# create_input_and_run_SCENT(37973,run_bs=TRUE,bootstrap_sig=TRUE)

print(paste0("Done running. Now processing ",nrow(chunkinfo)," peak-gene connections..."))

res.df.all = as.data.frame(do.call(rbind,out.sub))
res.df.all$fdr = NA
res.df.all$pval = res.df.all$boot_basic_p
res.df.all$pval[is.na(res.df.all$pval)] <- res.df.all$p[is.na(res.df.all$pval)]
res.df.all$fdr = p.adjust(res.df.all$pval,method = 'fdr')
res.df.all$celltype = celltype_to_use
res.df.all = res.df.all[order(res.df.all$i),]
res.df.all$i = 500*(num-1) + res.df.all$i
fwrite(res.df.all,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

print(paste0("Results written to: ",output_file))

