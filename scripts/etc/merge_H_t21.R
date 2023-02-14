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

numThreads = detectCores()/2
viridis = c("#FDE725FF","#440154FF")

########################################3
# Read in data: 

# Healthy data
dfseurat = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.ChromVAR.rds")
# subset to run things faster
# Ncell = 5000
# dfseurat = dfseurat[,sample(1:ncol(dfseurat),Ncell,replace = F)]
dfseurat$disease = "H"

dfseurat2 = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.ChromVAR.rds")
# subset to run things faster
# Ncell = 5000
# dfseurat2 = dfseurat2[,sample(1:ncol(dfseurat2),Ncell,replace = F)]
dfseurat2$disease = "T21"

########################################3

# combine H and T21 data:
dfseurat = merge(dfseurat,dfseurat2)
dfseurat = dfseurat[,!(dfseurat@meta.data$subclust_v6 %in% c("Unknown","No markers"))]

# remove empty peaks if necessary
rowcount <- rowSums(dfseurat@assays$ATAC@counts > 0)
rowKeep = rowcount > 0
dfseurat.sub = dfseurat.sub[rowKeep,]
# dfseurat@assays$ATAC@counts = (dfseurat@assays$ATAC@counts[rowKeep, ])
# dfseurat@assays$ATAC@ranges = (dfseurat@assays$ATAC@ranges[rowKeep,])

rm(dfseurat2)
###############

print("RunTFIDF...")
dfseurat <- RunTFIDF(dfseurat)
print("FindTopFeatures...")
dfseurat <- FindTopFeatures(dfseurat,min.cutoff = 'q0',assay="ATAC")
print("RunSVD...")
dfseurat <- RunSVD(dfseurat)
print("RunHarmony...")
lambda_val=1; dfseurat <- RunHarmony(dfseurat,"dataset",reduction="lsi",assay.use="ATAC",project.dim=FALSE,lambda=lambda_val,dims.use=2:30)
print("RunUMAP...")
dfseurat <- RunUMAP(dfseurat, reduction = "harmony",dims=1:50)

###############

umap_mat = dfseurat@reductions$umap@cell.embeddings
mutualknn30 <- getmutualknn(dfseurat@reductions$harmony@cell.embeddings, 30)
dfseurat@meta.data$disease2 = as.numeric(dfseurat@meta.data$disease=="T21")
dfseurat@meta.data$knn_disease_score = unlist(mclapply(1:ncol(mutualknn30),function(i) mean(dfseurat@meta.data$disease2[mutualknn30[,i]==1]),mc.cores = numThreads))

################

saveRDS(mutualknn30,file="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/shared_knn.rds")
saveRDS(dfseurat,"/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds")

#########################

tmp = data.frame(names=rownames(dfseurat@meta.data),
                 cell=dfseurat@meta.data$cell,
                 dataset=dfseurat@meta.data$dataset,
                 disease=dfseurat@meta.data$disease,
                 subclust_v6=dfseurat@meta.data$subclust_v6,
                 UMAP_1=umap_mat[,1],
                 UMAP_2=umap_mat[,2],
                 knn=dfseurat@meta.data$knn_disease_score)
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=knn),size=rel(0.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_viridis(begin=1,end=0)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.knn.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()

p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=disease),size=rel(0.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values=viridis)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.disease.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()

p2 <- ggplot(data=subset(tmp,disease=="H"), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=disease),size=1, na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") +
  scale_color_manual(values=viridis[1])
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.disease.h.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()

p2 <- ggplot(data=subset(tmp,disease=="T21"), aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=disease),size=1, na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") +
  scale_color_manual(values=viridis[2])
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.disease.T21.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()


library(RColorBrewer)
set.seed(03191995)
n <- length(unique(tmp$subclust_v6))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, n)
tmpout = aggregate(tmp[,c("UMAP_1","UMAP_2")],by=list(tmp$subclust_v6),mean)
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=subclust_v6),size=1, na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") +
  geom_richtext(data = tmpout,aes(label=Group.1,x=UMAP_1,y=UMAP_2)) +
  scale_color_manual(values=col_to_use)
# geom_richtext(label=tmp$Group.1,x=tmp$UMAP_1,y=tmp$UMAP_2)
# scale_color_gradientn(colors = viridis) + 

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.clusters.pdf")
pdf(f.plot,width=10.2,height=7)
print(p2)
dev.off()


