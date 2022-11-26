# module load R/4.1.2

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

#######################################################################################

peak_gene_links <- function(dfseurat) {
  
  print("Generating peaks coord file...")
  
  peaks = granges(dfseurat@assays$ATAC)
  peaks$peakid = rownames(dfseurat@assays$ATAC)
  peaks.df = as.data.frame(peaks@ranges)
  peaks.bed = data.frame(
    chr=as.character(peaks@seqnames),
    start=(peaks.df$start),
    end=(peaks.df$end),
    peakid=rownames(dfseurat@assays$ATAC)
  )
  
  print("Generating gene coord file...")
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
  
  print("Writing genes & peaks coord files...")
  fwrite(annot.bed,"~/tmp/genes.bed",quote = F,na = "NA",sep = "\t",row.names = F,col.names = F)
  fwrite(peaks.bed,"~/tmp/peaks.bed",quote = F,na = "NA",sep = "\t",row.names = F,col.names = F)
  
  print("Bedtools intersect genes & peaks coord files...")
  system("module load bedtools; bedtools intersect -a ~/tmp/genes.bed -b ~/tmp/peaks.bed -wa -wb > ~/tmp/genes_peaks_overlap.bed")
}

#######################################################################################
#######################################################################################
#######################################################################################

# Data input:
dfseurat=readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds")
dfseurat@meta.data$celltype=paste0(gsub(" ","_",dfseurat@meta.data$subclust_v6),"_",dfseurat@meta.data$disease)
f.cell.out = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/cells.txt")
fwrite(data.frame(unique(dfseurat@meta.data$celltype)),f.cell.out,row.names = F,col.names = F,quote = F,na = "NA",sep = "\t")
peak_gene_links(dfseurat)

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
  celltype=dfseurat@meta.data$celltype
)
rownames(meta) <- NULL

# extensions:
# - create interaction terms
# - run in a single cell type

# celltype_to_use = "HSCs_T21"
for (celltype_to_use in unique(dfseurat@meta.data$celltype)) {
  
  print(celltype_to_use)
  
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
  f.rna_out = paste0(fDir,"/rna/",celltype_to_use,".rna.rds")
  f.meta_out = paste0(fDir,"/meta/meta.",celltype_to_use,".rds")
  f.chunkinfo = paste0(fDir,"/chunkinfo/",celltype_to_use,".chunkinfo.txt")
  saveRDS(atac.all[peaks.keep,],f.atac_out)
  saveRDS(mrna[genes.keep,],f.rna_out)
  saveRDS(meta[cells.keep,],f.meta_out)
  fwrite(chunkinfo,f.chunkinfo,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
}


