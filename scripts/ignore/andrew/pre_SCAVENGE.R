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


########################################################################################
########################################################################################
########################################################################################

# This next analysis for HSCs:

# subset
ind = dfseurat@meta.data$subclust_v6 == "HSCs"
dfseurat.sub = dfseurat[,ind]

# remove empty peaks if necessary
rowcount <- rowSums(dfseurat.sub@assays$ATAC@counts > 0)
rowKeep = rowcount > 0
dfseurat.sub = dfseurat.sub[rowKeep,]
# dfseurat.sub = subset(dfseurat,rowKeep)
# dfseurat.sub@assays$ATAC@counts = (dfseurat.sub@assays$ATAC@counts[rowKeep, ])
# dfseurat.sub@assays$ATAC@ranges = (dfseurat.sub@assays$ATAC@ranges[rowKeep,])

print("RunTFIDF...")
dfseurat.sub <- RunTFIDF(dfseurat.sub)
print("FindTopFeatures...")
# dfseurat.sub <- FindTopFeatures(dfseurat.sub,min.cutoff = 'q0',assay="ATAC")
dfseurat.sub <- FindTopFeatures(dfseurat.sub,min.cutoff = 200,assay="ATAC")
print("RunSVD...")
dfseurat.sub <- RunSVD(dfseurat.sub)
print("RunHarmony...")
lambda_val=1; dfseurat.sub <- RunHarmony(dfseurat.sub,"dataset",reduction="lsi",assay.use="ATAC",project.dim=FALSE,lambda=lambda_val,dims.use=2:30)
print("RunUMAP...")
dfseurat.sub <- RunUMAP(dfseurat.sub, reduction = "harmony",dims=1:50)

###############

umap_mat = dfseurat.sub@reductions$umap@cell.embeddings
mutualknn30 <- getmutualknn(dfseurat.sub@reductions$harmony@cell.embeddings, 30)
dfseurat.sub@meta.data$disease2 = as.numeric(dfseurat.sub@meta.data$disease=="T21")
dfseurat.sub@meta.data$knn_disease_score = unlist(mclapply(1:ncol(mutualknn30),function(i) mean(dfseurat.sub@meta.data$disease2[mutualknn30[,i]==1]),mc.cores = numThreads))

############

tmp = data.frame(names=rownames(dfseurat.sub@meta.data),
                 cell=dfseurat.sub@meta.data$cell,
                 dataset=dfseurat.sub@meta.data$dataset,
                 disease=dfseurat.sub@meta.data$disease,
                 subclust_v6=dfseurat.sub@meta.data$subclust_v6,
                 UMAP_1=umap_mat[,1],
                 UMAP_2=umap_mat[,2],
                 knn=dfseurat.sub@meta.data$knn_disease_score)
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=knn),size=rel(0.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_viridis(begin=1,end=0)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.umap.knn.pdf")
pdf(f.plot,width=4.5,height=3.5)
print(p2)
dev.off()

p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=disease),size=rel(1.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values=viridis)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.umap.disease.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()


f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.umap.features.pdf")
pdf(f.plot,width=13,height=3.5)
DefaultAssay(dfseurat.sub) <- "RNA"
print(FeaturePlot(dfseurat.sub,features=c("CD34","SPINK2","PROM1"),ncol=3))
DefaultAssay(dfseurat.sub) <- "ATAC"
dev.off()

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.HSC_features.pdf")
pdf(f.plot,width=13,height=3.5)
DefaultAssay(dfseurat) <- "RNA"
print(FeaturePlot(dfseurat,features=c("CD34","SPINK2","PROM1"),ncol=3))
DefaultAssay(dfseurat) <- "ATAC"
dev.off()


tmp.full = data.frame(names=rownames(dfseurat@meta.data),
                      cell=dfseurat@meta.data$cell,
                      dataset=dfseurat@meta.data$dataset,
                      disease=dfseurat@meta.data$disease,
                      subclust_v6=dfseurat@meta.data$subclust_v6,
                      UMAP_1=dfseurat@reductions$umap@cell.embeddings[,1],
                      UMAP_2=dfseurat@reductions$umap@cell.embeddings[,2],
                      knn=dfseurat@meta.data$knn_disease_score)
tmp.sub = data.frame(names=rownames(dfseurat.sub@meta.data),
                     cell=dfseurat.sub@meta.data$cell,
                     dataset=dfseurat.sub@meta.data$dataset,
                     disease=dfseurat.sub@meta.data$disease,
                     subclust_v6=dfseurat.sub@meta.data$subclust_v6,
                     UMAP_1=dfseurat.sub@reductions$umap@cell.embeddings[,1],
                     UMAP_2=dfseurat.sub@reductions$umap@cell.embeddings[,2],
                     knn=dfseurat.sub@meta.data$knn_disease_score)
tmp.mg =  merge(tmp.full,tmp.sub[,c("names","knn")],by="names",all.x=TRUE)

tmp.mg1=subset(tmp.mg,is.na(knn.y))
tmp.mg2=subset(tmp.mg,!is.na(knn.y))
p2 <- ggplot(data=tmp.mg, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(data=tmp.mg1,aes(x=UMAP_1,y=UMAP_2),col='grey',size=rel(0.5), na.rm = TRUE, alpha = 0.6) +
  geom_point(data=tmp.mg2[sample(1:nrow(tmp.mg2),nrow(tmp.mg2),replace = F),],aes(x=UMAP_1,y=UMAP_2,color=disease),size=rel(0.5), na.rm = TRUE) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values=viridis)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.full_umap.disease.pdf")
pdf(f.plot,width=4.5,height=3.5)
print(p2)
dev.off()

p2 <- ggplot(data=tmp.mg, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(data=tmp.mg1,aes(x=UMAP_1,y=UMAP_2),col='grey',size=rel(0.5), na.rm = TRUE, alpha = 0.6) +
  geom_point(data=tmp.mg2[sample(1:nrow(tmp.mg2),nrow(tmp.mg2),replace = F),],aes(x=UMAP_1,y=UMAP_2,color=knn.y),size=rel(0.5), na.rm = TRUE) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_viridis(begin=1,end=0)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.full_umap.knn.pdf")
pdf(f.plot,width=4.5,height=3.5)
print(p2)
dev.off()


cor(tmp.mg$knn.x,tmp.mg$knn.y,use='na.or.complete')
table(tmp.mg$knn.x>0.8,tmp.mg$knn.y>0.8)

####################

saveRDS(dfseurat.sub,"/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

#############

DefaultAssay(dfseurat.sub) <- "RNA"
print("NormalizeData + FindVariableFeatures + ScaleData + RunPCA...")
dfseurat.sub <- NormalizeData(dfseurat.sub,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% # change to all genes using , features=rownames(pbmc) if using heatmap
  RunPCA()

print("RunHarmony...")
lambda_val=1; dfseurat.sub <- RunHarmony(dfseurat.sub,"dataset",lambda=lambda_val,reduction="pca",reduction.save="rna_harmony")

print("RunUMAP...")
dfseurat.sub <- RunUMAP(dfseurat.sub, reduction = "rna_harmony", dims = 1:30,reduction.key="rna_umap")

tmp = data.frame(names=rownames(dfseurat.sub@meta.data),
                      cell=dfseurat.sub@meta.data$cell,
                      dataset=dfseurat.sub@meta.data$dataset,
                      disease=dfseurat.sub@meta.data$disease,
                      subclust_v6=dfseurat.sub@meta.data$subclust_v6,
                      UMAP_1=dfseurat.sub@reductions$umap@cell.embeddings[,1],
                      UMAP_2=dfseurat.sub@reductions$umap@cell.embeddings[,2],
                      knn=dfseurat.sub@meta.data$knn_disease_score)
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=disease),size=rel(1.5), na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") + 
  scale_color_manual(values=viridis)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.disease.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.features.pdf")
pdf(f.plot,width=13,height=3.5)
DefaultAssay(dfseurat.sub) <- "RNA"
print(FeaturePlot(dfseurat.sub,features=c("CD34","SPINK2","PROM1"),ncol=3))
DefaultAssay(dfseurat.sub) <- "ATAC"
dev.off()

dfseurat.sub[["ATAC"]] <- NULL
dfseurat.sub@reductions$harmony <- NULL
saveRDS(dfseurat.sub,"/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA.HSC_only.rds")

dfseurat.sub.rna = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA.HSC_only.rds")
dfseurat.sub.atac = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")


f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.features2.pdf")
pdf(f.plot,width=14,height=8)
DefaultAssay(dfseurat.sub.rna) <- "RNA"
print(FeaturePlot(dfseurat.sub.rna,features=c("GATA1","KLF1","MPO","ITGA2B","MKI67","TOP2A"),ncol=3))
# DefaultAssay(dfseurat.sub.rna) <- "ATAC"
dev.off()

kclust=3
set.seed(031995)
kmeans_results <- kmeans(Embeddings(dfseurat.sub.rna,reduction="rna_harmony"),kclust)
kmeans_results1 <- data.frame(cluster=kmeans_results$cluster)
tmp = data.frame(names=rownames(dfseurat.sub.rna@meta.data),
                 cell=dfseurat.sub.rna@meta.data$cell,
                 dataset=dfseurat.sub.rna@meta.data$dataset,
                 disease=dfseurat.sub.rna@meta.data$disease,
                 subclust_v6=dfseurat.sub.rna@meta.data$subclust_v6,
                 UMAP_1=dfseurat.sub.rna@reductions$umap@cell.embeddings[,1],
                 UMAP_2=dfseurat.sub.rna@reductions$umap@cell.embeddings[,2],
                 knn=dfseurat.sub.rna@meta.data$knn_disease_score,
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
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.rna_umap.kmeans.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()


kclust=2
set.seed(031995)
# kmeans_results <- kmeans(Embeddings(dfseurat.sub.atac,reduction="harmony"),kclust)
# kmeans_results1 <- data.frame(cluster=kmeans_results$cluster)
tmp = data.frame(names=rownames(dfseurat.sub.atac@meta.data),
                 cell=dfseurat.sub.atac@meta.data$cell,
                 dataset=dfseurat.sub.atac@meta.data$dataset,
                 disease=dfseurat.sub.atac@meta.data$disease,
                 subclust_v6=dfseurat.sub.atac@meta.data$subclust_v6,
                 UMAP_1=dfseurat.sub.atac@reductions$umap@cell.embeddings[,1],
                 UMAP_2=dfseurat.sub.atac@reductions$umap@cell.embeddings[,2],
                 knn=dfseurat.sub.atac@meta.data$knn_disease_score,
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
# f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.atac_umap.kmeans.pdf")
# pdf(f.plot,width=9,height=7)
# print(p2)
# dev.off()
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.atac_umap.rna_kmeans.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()

kmeans_results <- kmeans(Embeddings(dfseurat.sub.rna,reduction="rna_harmony"),kclust)
traitName="mono"
f = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/old/scavenge/",traitName,".txt")
trait_mat = fread(f,data.table = F,stringsAsFactors = F) 
trait_mat$cell_dataset = paste(trait_mat$cell,trait_mat$dataset)
dfseurat.sub.rna@meta.data$i = 1:nrow(dfseurat.sub.rna@meta.data)
meta = merge(dfseurat.sub.rna@meta.data,trait_mat[,c("cell_dataset","TRS")],all.x=TRUE)
meta = meta[order(meta$i),]
rownames(meta) <- rownames(dfseurat.sub.rna@meta.data)
dfseurat.sub.rna@meta.data <- meta
kclust=3
set.seed(031995)
tmp = data.frame(UMAP_1=dfseurat.sub.rna@reductions$umap@cell.embeddings[,1],
                 UMAP_2=dfseurat.sub.rna@reductions$umap@cell.embeddings[,2],
                 disease=dfseurat.sub.rna@meta.data$disease,
                 TRS=dfseurat.sub.rna@meta.data$TRS,
                 kmeans=kmeans_results$cluster)
aggregate(tmp$TRS,by=list(cluster=tmp$kmeans,disease=tmp$disease),mean,na.rm=T)
aggregate(tmp$TRS,by=list(cluster=tmp$kmeans,disease=tmp$disease),median,na.rm=T)

table(tmp$kmeans)
# TRS U-MAP:
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2, color=TRS)) + 
  geom_point(size=1, na.rm = TRUE) + #, alpha = 0.6) +
  scale_color_gradientn(colors = viridis) + 
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.",traitName,".pdf")
pdf(f.plot,width=7,height=7)
print(p2)
dev.off()
aggregate(tmp$TRS,by=list(cluster=tmp$kmeans,disease=tmp$disease),mean,na.rm=T)
aggregate(tmp$TRS,by=list(cluster=tmp$kmeans,disease=tmp$disease),median,na.rm=T)

# TRS boxplot:
g=ggplot(tmp,aes(x=as.factor(kmeans),y=TRS,fill=disease)) +
  theme_bw() +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  labs(x="Cluster",y=paste0("TRS (",traitName,")"))
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_boxplots/HSCs.boxplot.",traitName,".pdf")
pdf(f.plot,width = 10,height=6)
print(g)
dev.off()

wilcox.test(tmp$TRS[tmp$kmeans==2],tmp$TRS[tmp$kmeans!=2])
t.test(tmp$TRS[tmp$kmeans==2],tmp$TRS[tmp$kmeans!=2])

f.out = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/RNA_kmeans.txt"
fwrite(tmp,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)



