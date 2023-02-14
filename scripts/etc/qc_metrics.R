Idents(df) = "orig.ident"
f.out <- paste0(dir,"/output/data/",DATASET,"_v2","/QC_metrics.h.pdf")
pdf(f.out,width=10,height = 4)
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","nCount_ATAC","TSS.enrichment"), ncol = 5)
dev.off()

DATASET="DS_Multiome_h"
df = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.ChromVAR.rds")

f.out = paste0(dir,"/output/data/",DATASET,"_v2","/meta.h.txt")
fwrite(df@meta.data,file = f.out,quote = F,na = "NA",sep = "\t",row.names = T,col.names = T)

DATASET="DS_Multiome_ds"
df = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.ChromVAR.rds")

f.out = paste0(dir,"/output/data/",DATASET,"_v2","/meta.ds.txt")
fwrite(df@meta.data,file = f.out,quote = F,na = "NA",sep = "\t",row.names = T,col.names = T)

library(data.table)
type='ds'
for (type in c("ds","h")) { 
  f=paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/meta.",type,".txt")
  df=fread(f,data.table = F,stringsAsFactors = F)
  df = df[,c("nFeature_RNA", "nCount_RNA", "percent.mt","nCount_ATAC","TSS.enrichment")]
  
  library(ggplot2)
  
  f.out = paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/qc.",type,".1.pdf")
  pdf(f.out,width=3,height = 4)
  g <- ggplot(df,aes(x=1,y=nFeature_RNA)) + labs(y="# expressed genes") + geom_violin(fill="turquoise") + geom_boxplot(width=0.1) + ggpubr::theme_pubr() + theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()); print(g)
  dev.off()
  
  f.out = paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/qc.",type,".2.pdf")
  pdf(f.out,width=3,height = 4)
  g <- ggplot(df,aes(x=1,y=log(nCount_RNA))) + labs(y="log(# RNA UMIs)") + geom_violin(fill="turquoise1") + geom_boxplot(width=0.1) + ggpubr::theme_pubr() + theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()); print(g)
  dev.off()
  
  f.out = paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/qc.",type,".3.pdf")
  pdf(f.out,width=3,height = 4)
  g <- ggplot(df,aes(x=1,y=percent.mt)) + labs(y="% mitochondrial reads") + geom_violin(fill="turquoise2") + geom_boxplot(width=0.07) + ggpubr::theme_pubr() + theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()); print(g)
  dev.off()
  
  f.out = paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/qc.",type,".4.pdf")
  pdf(f.out,width=3,height = 4)
  g <- ggplot(df,aes(x=1,y=log(nCount_ATAC))) + labs(y="log(# ATAC fragments in peaks)") + geom_violin(fill="turquoise3") + geom_boxplot(width=0.1) + ggpubr::theme_pubr() + theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()); print(g)
  dev.off()
  
  f.out = paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/qc.",type,".5.pdf")
  pdf(f.out,width=3,height = 4)
  g <- ggplot(df,aes(x=1,y=TSS.enrichment)) + labs(y="TSS enrichment score") + geom_violin(fill="turquoise4") + geom_boxplot(width=0.1) + ggpubr::theme_pubr() + theme(axis.title.x = element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()); print(g)
  dev.off()
  
}

