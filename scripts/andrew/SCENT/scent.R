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
options(stringsAsFactors = F)

## define functions

peak_gene_links <- function(dfseurat) {
  peaks = granges(dfseurat@assays$ATAC)
  peaks$peakid = rownames(dfseurat@assays$ATAC)
  peaks.bed = data.frame(
    chr=as.character(peaks@seqnames),
    start=(start(peaks)),
    end=(end(peaks)),
    peakid=rownames(dfseurat@assays$ATAC)
  )
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                                'start_position', 'end_position'),
                 filters = 'hgnc_symbol',
                 values = rownames(dfseurat@assays$RNA),
                 mart = ensembl)
  annot$start = annot$start-500000
  annot$end = annot$end+500000
  annot$chromosome_name = paste0("chr",annot$chromosome_name)
  annot.bed = annot[,c("chromosome_name","start_position","end_position","hgnc_symbol")]
  
  fwrite(annot.bed,"~/tmp/genes.bed",quote = F,na = "NA",sep = "\t",row.names = F,col.names = F)
  fwrite(peaks.bed,"~/tmp/peaks.bed",quote = F,na = "NA",sep = "\t",row.names = F,col.names = F)
  system("module load bedtools; bedtools intersect -a ~/tmp/genes.bed -b ~/tmp/peaks.bed -wa -wb > ~/tmp/genes_peaks_overlap.bed")
}

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

create_input_data = function(i) {
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

create_input_and_run_SCENT <- function(i,run_bs=TRUE) {
  if (i%%100==0) {print(i)}
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

# Data input:
dfseurat=readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")
celltype_to_use = "HSCs_H"

# peak_gene_links(dfseurat)
genes_peaks_overlap = fread("~/tmp/genes_peaks_overlap.bed",data.table = F,stringsAsFactors = F)
genes_peaks_overlap = genes_peaks_overlap[,c("V4","V8")]

# meta data:
meta = data.frame(
  cell=rownames(dfseurat@meta.data),
  percent_mito = dfseurat@meta.data$percent.mt,
  nUMI=colSums(dfseurat@assays$RNA@counts),
  sample=dfseurat@meta.data$dataset,
  celltype0 = gsub(" ","_",dfseurat@meta.data$subclust_v6),
  disease0 = dfseurat@meta.data$disease,
  celltype=paste0(gsub(" ","_",dfseurat@meta.data$subclust_v6),"_",dfseurat@meta.data$disease)
)
rownames(meta) <- NULL

# extend to interaction model? with run in cell type

cells.keep = paste0(gsub(" ","_",dfseurat@meta.data$subclust_v6),"_",dfseurat@meta.data$disease) == celltype_to_use

pct.rna.keep = 0.05
pct.atac.keep = 0.05
atac.all <- dfseurat@assays$ATAC@counts[,cells.keep]
atac.all@x[atac.all@x > 1] <- 1

pct.accessible = rowSums(atac.all)/ncol(atac.all)
peak.data = data.frame(peakid=rownames(atac.all),pct=pct.accessible)
ind.peaks.keep = pct.accessible > pct.atac.keep
peaks.keep = rownames(atac.all)[ind.peaks.keep]
# atac = atac.all[peaks.keep,]

mrna <- dfseurat@assays$RNA@counts[,cells.keep]
mrna.binarized = mrna
mrna.binarized@x[mrna.binarized@x > 1] <- 1
pct.expressed = rowSums(mrna.binarized)/ncol(mrna.binarized)
ind.genes.keep = pct.expressed > pct.rna.keep
genes.keep = rownames(mrna)[ind.genes.keep]

# extract genes and peaks that are to be tested:
chunkinfo = genes_peaks_overlap
colnames(chunkinfo)<-c("gene","peak")
chunkinfo = subset(chunkinfo,gene %in% genes.keep & peak %in% peaks.keep)

fDir = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input"
f.atac_out = paste0(fDir,"/atac/",celltype_to_use,".atac.rds")
f.rna_out = paste0(fDir,"/",celltype_to_use,".rna.rds")
f.meta_out = paste0(fDir,"/meta.rds")
f.chunkinfo = paste0(fDir,"/",celltype_to_use,".chunkinfo.txt")
saveRDS(atac.all[peaks.keep,],f.atac_out)
saveRDS(mrna[genes.keep,],f.rna_out)
saveRDS(meta,f.meta_out)
fwrite(chunkinfo,f.chunkinfo,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

######################

atac.all = readRDS(f.atac_out)
mrna = readRDS(f.rna_out)
meta = readRDS(f.meta_out)
chunkinfo = fread(chunkinfo,data.table = F,stringsAsFactors = F)

numThreads = detectCores()/2
out = mclapply(1:nrow(chunkinfo),create_input_and_run_SCENT,run_bs=FALSE,mc.cores = numThreads)
res.df = as.data.frame(do.call(rbind,out))
res.df$fdr = p.adjust(res.df$p,method='fdr')
output_file = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
fwrite(res.df,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
res.df[order(res.df$fdr)[1:4],]; mean(res.df$fdr < 0.1)
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



###############

# celltype_to_use = "HSCs_H"
# res<-list()
# for (i in 12:nrow(chunkinfo)) {
#   print(i)
#   # i = 1
#   # extract data and create dataframe:
#   
#   if(!(nonzero_m > pct.rna.keep & nonzero_a > pct.atac.keep)) { print(paste0(i," error!"))}
#   
#   if(nonzero_m > pct.rna.keep & nonzero_a > pct.atac.keep){
#     # poisson
#     base = glm(exprs ~ ., family = 'poisson', data = df2.input) # why do raw counts instead of normalization? or voom re-weighting?
#     coefs<-summary(base)$coefficients["atac",]
#     bs = boot::boot(df2.input,assoc_poisson, R = 100, stype = 'i')
#     p0 = basic_p(bs$t0[1], bs$t[,1])
#     if(p0<0.1){
#       bs = boot::boot(df2.input,assoc_poisson, R = 500, stype = 'i')
#       p0 = basic_p(bs$t0[1], bs$t[,1])
#     }
#     if(p0<0.05){
#       bs = boot::boot(df2.input,assoc_poisson, R = 2500, stype = 'i')
#       p0 = basic_p(bs$t0[1], bs$t[,1])
#     }
#     if(p0<0.01){
#       bs = boot::boot(df2.input,assoc_poisson, R = 25000, stype = 'i')
#       p0 = basic_p(bs$t0[1], bs$t[,1])
#     }
#     if(p0<0.001){
#       bs = boot::boot(df2.input,assoc_poisson, R = 50000, stype = 'i')
#       p0 = basic_p(bs$t0[1], bs$t[,1])
#     }
#     out <- data.frame(i=i,gene=gene,peak=this_peak,beta=coefs[1],se=coefs[2],z=coefs[3],p=coefs[4],boot_basic_p=p0)
#     res[[i]]<-out
#   }
# }
# res.df = as.data.frame(do.call(rbind,res))
# res.df$celltype = celltype_to_use
# fwrite(res.df,output_file,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
# write.table(res, output, quote=F, row=F, sep="\t")






