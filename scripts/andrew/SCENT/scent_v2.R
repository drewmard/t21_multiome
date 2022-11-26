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
f.cell.out = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/cells.txt")
cells = fread(f.cell.out,data.table = F,stringsAsFactors = F,header = F)

for (celltype_to_use in cells[,1]) {

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

numThreads = detectCores() #/2

print(paste0("Running ",nrow(chunkinfo)," peak-gene connections..."))

# First pass poisson regression model:
print(Sys.time())
out = mclapply(1:nrow(chunkinfo),create_input_and_run_SCENT,run_bs=FALSE,mc.cores = numThreads)
print(Sys.time())
res.df = as.data.frame(do.call(rbind,out))
res.df$fdr = p.adjust(res.df$p,method='fdr')
output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
fwrite(res.df,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
print(res.df[order(res.df$fdr)[1:4],]); print(mean(res.df$fdr < 0.1))
}

# Second pass poisson regression model:
numCores_to_use = min(sum(res.df$fdr<0.1)+1,numThreads)
out.sub = mclapply(which(res.df$fdr<0.1),create_input_and_run_SCENT,
                   run_bs=TRUE,
                   bootstrapping_sig=TRUE,
                   mc.cores = numCores_to_use)
res.df.sub = as.data.frame(do.call(rbind,out.sub))
res.df.sub$fdr = NA
res.df.all = subset(as.data.frame(rbind(res.df.sub,res.df)),!duplicated(i))
res.df.all$pval = res.df.all$boot_basic_p
res.df.all$pval[is.na(res.df.all$pval)] <- res.df.all$p[is.na(res.df.all$pval)]
res.df.all$fdr = p.adjust(res.df.all$pval,method = 'fdr')
res.df.all$celltype = celltype_to_use
output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out2/",celltype_to_use,".txt")
fwrite(res.df,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
