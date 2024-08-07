# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# module load R/4.1.2

# Code adapted from SCENT: https://github.com/immunogenomics/SCENT/blob/main/SCENT.R (original author: Saori Sakaue)

library(biomaRt)
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

#######################################################################################

## define functions

#' Interpolate a p-value from quantiles that should be "null scaled"
#'
#' @param q bootstrap quantiles, centered so that under the null, theta = 0
#' @return two-sided p-value
interp_pval = function(q) {
  R = length(q)
  tstar = sort(q)
  zero = findInterval(0, tstar)
  if(zero == 0 || zero == R) return(2/R) # at/beyond extreme values
  pval = 2*min(zero/R, (R-zero)/R)
  pval
}

#' Derive a p-value from a vector of bootstrap samples using the "basic" calculation
#'
#' @param obs observed value of parameter (using actual data)
#' @param boot vector of bootstraps
#'
#' @return p-value
basic_p = function(obs, boot, null = 0){
  interp_pval(2*obs - boot - null)
}


assoc_poisson = function(data, idx = seq_len(nrow(data))){
  gg = glm(exprs ~ ., family = 'poisson', data = data[idx,,drop = FALSE])
  # gg = glm(exprs ~ atac + percent_mito + log(nUMI) + sample, family = 'poisson', data = data[idx,,drop = FALSE])
  c(coef(gg)['atac'], diag(vcov(gg))['atac'])
}

create_input_data = function(i,pct.rna.keep=0.05,pct.atac.keep=0.05) {
  gene=chunkinfo$gene[i]
  this_peak=chunkinfo$peak[i]
  atac_target<-data.frame(cell=colnames(atac.all),atac=as.numeric(atac.all[this_peak,]))
  mrna_target<-mrna[gene,]
  df <- data.frame(cell=names(mrna_target),exprs=as.numeric(mrna_target))
  df<-merge(df,atac_target,by="cell")
  df<-merge(df,meta,by="cell")
  
  # Subset cells to test:
  df2 <- df[df$celltype==celltype_to_use,]
  
  # Binarize peaks:
  df2[df2$atac>0,]$atac<-1
  
  # QC: Require >5% expressed genes and accessible peaks:
  nonzero_m  <- length( df2$exprs[ df2$exprs > 0] ) / length( df2$exprs )
  nonzero_a  <- length( df2$atac[ df2$atac > 0] ) / length( df2$atac )
  if(!(nonzero_m > pct.rna.keep & nonzero_a > pct.atac.keep)){return(NULL)}
  # create log nUMI column
  df2$log_nUMI = log(df2$nUMI)
  
  df2.input = df2[,c("exprs","atac","percent_mito","log_nUMI","sample")]
  return(df2.input)
}

bootstrapping = function(df2.input,i) {
  bs = boot::boot(df2.input,assoc_poisson, R = 100, stype = 'i')
  p0 = basic_p(bs$t0[1], bs$t[,1])
  print(paste0("R=100, i=",i,": p0=",p0))
  if(p0<0.1){
    bs = boot::boot(df2.input,assoc_poisson, R = 500, stype = 'i')
    p0 = basic_p(bs$t0[1], bs$t[,1])
    print(paste0("R=500, i=",i,": p0=",p0))
  }
  if(p0<0.05){
    bs = boot::boot(df2.input,assoc_poisson, R = 2500, stype = 'i')
    p0 = basic_p(bs$t0[1], bs$t[,1])
    print(paste0("R=2500, i=",i,": p0=",p0))
  }
  if(p0<0.01){
    bs = boot::boot(df2.input,assoc_poisson, R = 25000, stype = 'i')
    p0 = basic_p(bs$t0[1], bs$t[,1])
    print(paste0("R=25000, i=",i,": p0=",p0))
  }
  if(p0<0.001){
    bs = boot::boot(df2.input,assoc_poisson, R = 50000, stype = 'i')
    p0 = basic_p(bs$t0[1], bs$t[,1])
    print(paste0("R=50000, i=",i,": p0=",p0))
  }
  return(p0)
}

bootstrapping_sig = function(df2.input,i) {
  bs = boot::boot(df2.input,assoc_poisson, R = 50000, stype = 'i')
  p0 = basic_p(bs$t0[1], bs$t[,1])
  print(paste0("R=50000, i=",i,": p0=",p0))
  return(p0)
}

SCENT = function(df2.input,i,run_bs=TRUE,bootstrap_sig=FALSE) {
  # poisson
  base = glm(exprs ~ ., family = 'poisson', data = df2.input) # why do raw counts instead of normalization? or voom re-weighting?
  coefs<-summary(base)$coefficients["atac",]
  if (run_bs) {
    if (bootstrap_sig) {
      p0 = bootstrapping_sig(df2.input,i)
    } else {
      p0 = bootstrapping(df2.input,i)
    }
  } else {
    p0 = NA
  }
  # return data:
  gene=chunkinfo$gene[i]
  this_peak=chunkinfo$peak[i]
  out <- data.frame(i=i,gene=gene,peak=this_peak,beta=coefs[1],se=coefs[2],z=coefs[3],p=coefs[4],boot_basic_p=p0)
  return(out)
}

create_input_and_run_SCENT <- function(i,run_bs=TRUE,iter_print=1000) {
  if ( (i %% iter_print) == 0) {print(i)} # print update i every $iter_print iterations
  df2.input = create_input_data(i = i)
  if(!is.null(df2.input)){
    out = SCENT(df2.input = df2.input,i=i,run_bs=run_bs)
    return(out)
  } else {
    return(data.frame())
  }
}

###############################################################################
###############################################################################
###############################################################################

# celltype_to_use = "HSCs_T21"
library(data.table)

peak_of_interest="chr6-135096824-135097623"
gene_of_interest = "HBS1L"

celltype_to_use = "HSCs_H"
print(celltype_to_use)
fDir = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input"
f.atac_out = paste0(fDir,"/atac/",celltype_to_use,".atac.rds")
f.rna_out = paste0(fDir,"/rna/",celltype_to_use,".rna.rds")
f.meta_out = paste0(fDir,"/meta/meta.",celltype_to_use,".rds")
f.chunkinfo = paste0(fDir,"/chunkinfo/",celltype_to_use,".chunkinfo.txt")

atac.all = readRDS(f.atac_out)
mrna = readRDS(f.rna_out)
meta = readRDS(f.meta_out)
chunkinfo = fread(f.chunkinfo,data.table = F,stringsAsFactors = F)

gene=gene_of_interest
this_peak=peak_of_interest
dfseurat=readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

atac_target<-data.frame(cell=colnames(atac.all),atac=as.numeric(atac.all[this_peak,]))
mrna_target<-mrna[gene,]
df <- data.frame(cell=names(mrna_target),exprs=as.numeric(mrna_target))
df<-merge(df,atac_target,by="cell")
df<-merge(df,meta,by="cell")

# Subset cells to test:
df2 <- df[df$celltype==celltype_to_use,]

# Binarize peaks:
df2[df2$atac>0,]$atac<-1

# QC: Require >5% expressed genes and accessible peaks:
nonzero_m  <- length( df2$exprs[ df2$exprs > 0] ) / length( df2$exprs )
nonzero_a  <- length( df2$atac[ df2$atac > 0] ) / length( df2$atac )

