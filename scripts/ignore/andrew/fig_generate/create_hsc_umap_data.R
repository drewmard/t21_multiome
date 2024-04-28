# Run on cluster!
# This file creates a UMAP file needed for figure panels L-P.

library(Seurat)
library(Signac)

dfseurat.sub.rna = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA.HSC_only.rds")
dfseurat.sub.atac = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

# which(rownames(dfseurat.sub.atac@meta.data) != rownames(dfseurat.sub.rna@meta.data))

tmp = data.frame(names=rownames(dfseurat.sub.rna@meta.data),
                 cell=dfseurat.sub.rna@meta.data$cell,
                 dataset=dfseurat.sub.rna@meta.data$dataset,
                 disease=dfseurat.sub.rna@meta.data$disease,
                 subclust_v6=dfseurat.sub.rna@meta.data$subclust_v6,
                 UMAP_1_RNA=dfseurat.sub.rna@reductions$umap@cell.embeddings[,1],
                 UMAP_2_RNA=dfseurat.sub.rna@reductions$umap@cell.embeddings[,2],
                 UMAP_1_ATAC=dfseurat.sub.atac@reductions$umap@cell.embeddings[,1],
                 UMAP_2_ATAC=dfseurat.sub.atac@reductions$umap@cell.embeddings[,2],
                 knn_ATAC=dfseurat.sub.atac@meta.data$knn_disease_score,
                 kmeans_RNA=kmeans_results$cluster)

# aggregate(knn_ATAC ~ as.factor(kmeans_RNA),tmp,mean)
# aggregate(disease=="T21" ~ as.factor(kmeans_RNA),tmp,mean)

fwrite(tmp,"/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/hsc.kmeans.UMAP.txt",row.names = F,col.names = T,quote = F,na = "NA",sep = "\t")
