library(Seurat)
library(Signac)
library(data.table)

# Data input:
dfseurat=readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds")
dfseurat@meta.data$celltype=paste0(gsub(" ","_",dfseurat@meta.data$subclust_v6),"_",dfseurat@meta.data$disease)

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

celltype_to_use = "HSCs_T21"
for (celltype_to_use in c("HSCs_H","HSCs_T21")) {
  
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
expr.data = data.frame(geneid=rownames(mrna.binarized),pct=pct.expressed)

fout = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out_info/expr.",celltype_to_use,".txt")
fwrite(expr.data,fout,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
fout = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out_info/peak.",celltype_to_use,".txt")
fwrite(peak.data,fout,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

