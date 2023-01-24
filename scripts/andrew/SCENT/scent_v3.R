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

downsample=FALSE

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

bootstrapping_sig = function(df2.input,i,pval,res.df=NULL) {
  if (pval >= 0.1) {
    bs = boot::boot(df2.input,assoc_poisson, R = 100, stype = 'i')
    p0 = basic_p(bs$t0[1], bs$t[,1])
    if (p0 < 0.1) {
      bs = boot::boot(df2.input,assoc_poisson, R = 1000, stype = 'i')
      p0 = basic_p(bs$t0[1], bs$t[,1])
    }
    if (p0 < 0.01) {
      bs = boot::boot(df2.input,assoc_poisson, R = 10000, stype = 'i')
      p0 = basic_p(bs$t0[1], bs$t[,1])
    }
  } else if (pval >= 0.01 & pval < 0.1) {
    bs = boot::boot(df2.input,assoc_poisson, R = 1000, stype = 'i')
    p0 = basic_p(bs$t0[1], bs$t[,1])
    if (p0 < 0.01) {
      bs = boot::boot(df2.input,assoc_poisson, R = 10000, stype = 'i')
      p0 = basic_p(bs$t0[1], bs$t[,1])
    }
  } else {
    bs = boot::boot(df2.input,assoc_poisson, R = 10000, stype = 'i')
    p0 = basic_p(bs$t0[1], bs$t[,1])
  }
  return(p0)
}

SCENT = function(df2.input,i,run_bs=TRUE,bootstrap_sig=FALSE) {
  # poisson
  base = glm(exprs ~ ., family = 'poisson', data = df2.input) # why do raw counts instead of normalization? or voom re-weighting?
  coefs<-summary(base)$coefficients["atac",]
  if (run_bs) {
    if (bootstrap_sig) {
      p0 = bootstrapping_sig(df2.input,i,pval=coefs[4])
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

create_input_and_run_SCENT <- function(i,run_bs=TRUE,bootstrap_sig=FALSE,iter_print=1000) {
  if (run_bs) {print(i)}
  if ( (i %% iter_print) == 0) {print(i)} # print update i every $iter_print iterations
  df2.input = create_input_data(i = i)
  if(!is.null(df2.input)){
    out = SCENT(df2.input = df2.input,i=i,run_bs=run_bs,bootstrap_sig=bootstrap_sig)
    return(out)
  } else {
    return(data.frame())
  }
}

###############################################################################
###############################################################################
###############################################################################

# celltype_to_use = "HSCs_T21"
# f.cell.out = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/cells.txt")
# cells = fread(f.cell.out,data.table = F,stringsAsFactors = F,header = F)[,1]
# cells = cells[c(7,21,(1:length(cells))[-c(7,21)])]

# celltype_to_use = "HSCs_H"
# for (celltype_to_use in cells) {

# source("/oak/stanford/groups/smontgom/amarder/t21_multiome/scripts/andrew/SCENT/scent_v4")

projName="tmparm"
celltype_to_use = "HSCs_T21"
num=77

args = commandArgs(trailingOnly=TRUE)
projName = args[1]
celltype_to_use = args[2]
num = as.numeric(args[3])

print(celltype_to_use)
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

print(paste0("Running ",nrow(chunkinfo)," peak-gene connections..."))
system.time(out.sub <- lapply(1:3, #nrow(chunkinfo),
                              create_input_and_run_SCENT,
                              run_bs=TRUE,
                              bootstrap_sig=TRUE))
numThreads = detectCores() #/2
# system.time(out.sub <- mclapply(1:3,
#                                 create_input_and_run_SCENT,
#                                 run_bs=TRUE,
#                                 bootstrap_sig=TRUE,
#                                 mc.cores = min(nrow(chunkinfo),numThreads)
#                                 ))
# create_input_and_run_SCENT(37973,run_bs=TRUE,bootstrap_sig=TRUE)
res.df.all = as.data.frame(do.call(rbind,out.sub))
res.df.all$fdr = NA
res.df.all$pval = res.df.all$boot_basic_p
res.df.all$pval[is.na(res.df.all$pval)] <- res.df.all$p[is.na(res.df.all$pval)]
res.df.all$fdr = p.adjust(res.df.all$pval,method = 'fdr')
res.df.all$celltype = celltype_to_use
res.df.all = res.df.all[order(res.df.all$i),]
res.df.all$i = 500*(num-1) + res.df.all$i
fwrite(res.df.all,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)










if (downsample) {
  set.seed(03191995)
  ind = sort(sample(1:3784,2431,replace = F))
  atac.all = atac.all[,ind]
  mrna = mrna[,ind]
  meta = meta[ind,]
  chunkinfo = fread(f.chunkinfo,data.table = F,stringsAsFactors = F)
  
  # First pass poisson regression model:
  print(Sys.time())
  out = mclapply(1:nrow(chunkinfo),create_input_and_run_SCENT,run_bs=FALSE,mc.cores = numThreads)
  print(Sys.time())
  res.df = as.data.frame(do.call(rbind,out))
  res.df$fdr = p.adjust(res.df$p,method='fdr')
  output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out1/",celltype_to_use,".down.txt")
  fwrite(res.df,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  print(res.df[order(res.df$fdr)[1:4],]); print(mean(res.df$fdr < 0.1))
  
  numThreads = detectCores() #/2
  
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out1/",celltype_to_use,".down.txt")
  res.df = fread(f,data.table = F,stringsAsFactors = F)
  
  # Second pass poisson regression model:
  numCores_to_use = numThreads # min(sum(res.df$fdr<0.1)+1,numThreads)
  # out.sub = mclapply(which(res.df$fdr<0.1),create_input_and_run_SCENT,
  #                    run_bs=TRUE,
  #                    bootstrap_sig=TRUE,
  #                    res.df.input = res.df,
  #                    mc.cores = numCores_to_use)
  out.sub = mclapply(1:nrow(chunkinfo),create_input_and_run_SCENT,
                     run_bs=TRUE,
                     bootstrap_sig=TRUE,
                     res.df.input = res.df,
                     mc.cores = numCores_to_use)
  out.sub = mclapply(1:20,create_input_and_run_SCENT,
                     run_bs=TRUE,
                     bootstrap_sig=TRUE,
                     res.df.input = res.df,
                     mc.cores = numCores_to_use)
  # create_input_and_run_SCENT(37973,run_bs=TRUE,bootstrap_sig=TRUE)
  res.df.sub = as.data.frame(do.call(rbind,out.sub))
  res.df.sub$fdr = NA
  res.df.all = subset(as.data.frame(rbind(res.df.sub,res.df)),!duplicated(i))
  res.df.all$pval = res.df.all$boot_basic_p
  res.df.all$pval[is.na(res.df.all$pval)] <- res.df.all$p[is.na(res.df.all$pval)]
  res.df.all$fdr = p.adjust(res.df.all$pval,method = 'fdr')
  res.df.all$celltype = celltype_to_use
  res.df.all = res.df.all[order(res.df.all$i),]
  output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out2/",celltype_to_use,".down.txt")
  fwrite(res.df.all,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
} else {
  
  # First pass poisson regression model:
  print(Sys.time())
  out = mclapply(1:nrow(chunkinfo),create_input_and_run_SCENT,run_bs=FALSE,mc.cores = numThreads)
  print(Sys.time())
  res.df = as.data.frame(do.call(rbind,out))
  res.df$fdr = p.adjust(res.df$p,method='fdr')
  output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
  fwrite(res.df,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  print(res.df[order(res.df$fdr)[1:4],]); print(mean(res.df$fdr < 0.1))
  
  
  library(data.table)
  library(parallel)
  numThreads = detectCores() #/2
  
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
  res.df = fread(f,data.table = F,stringsAsFactors = F)
  
  # Second pass poisson regression model:
  numCores_to_use = min(sum(res.df$fdr<0.2)+1,numThreads)
  out.sub = mclapply(which(res.df$fdr<0.2),create_input_and_run_SCENT,
                     run_bs=TRUE,
                     bootstrap_sig=TRUE,
                     mc.cores = numCores_to_use)
  # create_input_and_run_SCENT(37973,run_bs=TRUE,bootstrap_sig=TRUE)
  res.df.sub = as.data.frame(do.call(rbind,out.sub))
  res.df.sub$fdr = NA
  res.df.all = subset(as.data.frame(rbind(res.df.sub,res.df)),!duplicated(i))
  res.df.all$pval = res.df.all$boot_basic_p
  res.df.all$pval[is.na(res.df.all$pval)] <- res.df.all$p[is.na(res.df.all$pval)]
  res.df.all$fdr = p.adjust(res.df.all$pval,method = 'fdr')
  res.df.all$celltype = celltype_to_use
  res.df.all = res.df.all[order(res.df.all$i),]
  output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out2/",celltype_to_use,".txt")
  fwrite(res.df.all,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
}