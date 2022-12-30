library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(gchromVAR)
library(SCAVENGE)
library(harmony)
library(dplyr)
library(uwot)
library(parallel)
library(ggplot2)
library(data.table)
library(viridis)
library(ggtext)

dfseurat.sub.rna = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA.HSC_only.rds")
kclust=3
set.seed(031995)
kmeans_results <- kmeans(Embeddings(dfseurat.sub.rna,reduction="rna_harmony"),kclust)
dfseurat.sub.rna@meta.data$kmeans = kmeans_results$cluster

Idents(dfseurat.sub.rna) <- "kmeans"
de.out <- FindAllMarkers(dfseurat.sub.rna,logfc.threshold=1,min.pct = 0.2)
de.out

ind = grep("^MT-",rownames(dfseurat.sub.rna))
rownames(dfseurat.sub.rna)[ind]
dfseurat.sub.rna = dfseurat.sub.rna[-ind,]

dfseurat.sub.rna@reductions$pca <- NULL
dfseurat.sub.rna@reductions$rna_harmony <- NULL
dfseurat.sub.rna@reductions$umap <- NULL
DefaultAssay(dfseurat.sub.rna) <- "RNA"
print("NormalizeData + FindVariableFeatures + ScaleData + RunPCA...")
dfseurat.sub.rna <- NormalizeData(dfseurat.sub.rna,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% # change to all genes using , features=rownames(pbmc) if using heatmap
  RunPCA()

print("RunHarmony...")
lambda_val=1; dfseurat.sub.rna <- RunHarmony(dfseurat.sub.rna,"dataset",lambda=lambda_val,reduction="pca",reduction.save="harmony")

print("RunUMAP...")
dfseurat.sub.rna <- RunUMAP(dfseurat.sub.rna, reduction = "harmony", dims = 1:30,reduction.key="umap")

kclust=3
set.seed(031995)
kmeans_results <- kmeans(Embeddings(dfseurat.sub.rna,reduction="harmony"),kclust)
dfseurat.sub.rna@meta.data$kmeans = kmeans_results$cluster

Idents(dfseurat.sub.rna) <- "kmeans"
de.out <- FindAllMarkers(dfseurat.sub.rna,logfc.threshold=2,min.pct = 0.2)
de.out

table(kmeans_results$cluster,dfseurat.sub.rna@meta.data$disease)

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.no_mito.features.pdf")
pdf(f.plot,width=13,height=7)
DefaultAssay(dfseurat.sub.rna) <- "RNA"
print(FeaturePlot(dfseurat.sub.rna,features=c("ANK1","SLC25A21","IL7","ASPM","BRIP1","ALB"),ncol=3))
# DefaultAssay(dfseurat.sub) <- "ATAC"
dev.off()

tmp = data.frame(names=rownames(dfseurat.sub.rna@meta.data),
                 cell=dfseurat.sub.rna@meta.data$cell,
                 dataset=dfseurat.sub.rna@meta.data$dataset,
                 disease=dfseurat.sub.rna@meta.data$disease,
                 subclust_v6=dfseurat.sub.rna@meta.data$subclust_v6,
                 UMAP_1=dfseurat.sub.rna@reductions$umap@cell.embeddings[,1],
                 UMAP_2=dfseurat.sub.rna@reductions$umap@cell.embeddings[,2],
                 kmeans=kmeans_results$cluster)
library(RColorBrewer)
set.seed(031995)
n <- length(unique(tmp$kmeans))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, n)
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=as.factor(kmeans)),size=rel(1.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + labs(col="cluster") +
  scale_color_manual(values=col_to_use)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.no_mito.kmeans.pdf")
pdf(f.plot,width=9*0.7,height=7*0.7)
print(p2)
dev.off()

viridis = c("#FDE725FF","#440154FF")
viridis = c("yellow3","#440154FF")

p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=disease),size=rel(0.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values=viridis)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.no_mito.disease.pdf")
pdf(f.plot,width=9*0.7,height=7*0.7)
print(p2)
dev.off()

f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/FindAllMarkers.kmeans.txt")
print(paste0("Saving to file: ",f.out))
print("...")
fwrite(de.out,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.no_mito.qc_features.pdf")
pdf(f.plot,width=8,height=3.5)
DefaultAssay(dfseurat.sub.rna) <- "RNA"
print(FeaturePlot(dfseurat.sub.rna,features=c("nCount_RNA","nFeature_RNA"),ncol=2))
# DefaultAssay(dfseurat.sub) <- "ATAC"
dev.off()

aggregate(dfseurat.sub.rna@meta.data[,c("nCount_RNA","nFeature_RNA")],by=list(disease=dfseurat.sub.rna@meta.data$disease,kmeans=dfseurat.sub.rna@meta.data$kmeans),median)




