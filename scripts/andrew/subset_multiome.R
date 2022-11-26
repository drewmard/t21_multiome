# subset_multiome

dfrna = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/Multiome.RNA_ATAC.h.rds")
dfrna$disease = "H"
DefaultAssay(dfrna) <- "RNA"
dfrna[["ATAC"]] <- NULL

dfrna2 = readRDS("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.ds.rds")
dfrna2$disease = "T21"
DefaultAssay(dfrna2) <- "RNA"
dfrna2[["ATAC"]] <- NULL

dfrna = merge(dfrna,dfrna2)
dfrna = dfrna[,(dfrna$cell_dataset %in% dfseurat$cell_dataset)]
which(rownames(dfrna@meta.data)!=rownames(dfseurat@meta.data))
dfseurat[["RNA"]] <- dfrna[["RNA"]]

dfrna.sub = dfrna[,(dfrna$cell_dataset %in% dfseurat.sub$cell_dataset)]
which(rownames(dfrna.sub@meta.data)!=rownames(dfseurat.sub@meta.data))
dfseurat.sub[["RNA"]] <- dfrna.sub[["RNA"]]

dfseurat.sub@assays$RNA@counts[1:10,1:10]

######

DefaultAssay(dfseurat) <- "RNA"
dfseurat[["ATAC"]] <- NULL
# dfseurat[["chromvar"]] <- NULL
ind = dfseurat@meta.data$subclust_v6 %in% c("HSCs","MEMPs")
dfseurat.sub = dfseurat[,ind]
meta = dfseurat.sub@meta.data
meta$i = 1:nrow(meta)
meta = merge(meta,tmp[,c("UMAP_1","kmeans")],by="row.names",all.x=TRUE)
meta = meta[order(meta$i),]
meta$kmeans[is.na(meta$kmeans)] <- 4
dfseurat.sub@meta.data <- meta

DefaultAssay(dfseurat.sub) <- "RNA"
print("NormalizeData + FindVariableFeatures + ScaleData + RunPCA...")
dfseurat.sub <- NormalizeData(dfseurat.sub,normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% 
  ScaleData() %>% # change to all genes using , features=rownames(pbmc) if using heatmap
  RunPCA()

print("RunHarmony...")
# lambda_val=1; dfseurat.sub <- RunHarmony(dfseurat.sub,"dataset",lambda=lambda_val,reduction="pca",reduction.save="_rna_harmony")
dfseurat.sub <- RunHarmony(dfseurat.sub,"dataset",lambda=1,reduction="pca")
dfseurat.sub@reductions$rna_harmony <- dfseurat.sub@reductions$harmony

print("RunUMAP...")
dfseurat.sub <- RunUMAP(dfseurat.sub, reduction = "rna_harmony", dims = 1:30)

tmp = data.frame(names=rownames(dfseurat.sub@meta.data),
                 cell=dfseurat.sub@meta.data$cell,
                 dataset=dfseurat.sub@meta.data$dataset,
                 disease=dfseurat.sub@meta.data$disease,
                 subclust_v6=dfseurat.sub@meta.data$subclust_v6,
                 UMAP_1=dfseurat.sub@reductions$umap@cell.embeddings[,1],
                 UMAP_2=dfseurat.sub@reductions$umap@cell.embeddings[,2],
                 kmeans=dfseurat.sub@meta.data$kmeans)
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
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs_MEMPs.rna_umap.kmeans.pdf")
pdf(f.plot,width=9,height=7)
print(p2)
dev.off()

dir.create("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/violin")
Idents(dfseurat.sub) <- "kmeans"
genes.sub = c("PROM1", "SPINK2", "GATA1", "ITGA2B" , "MKI67","PROCR")
rng = seq(1,length(genes.sub));
for (k in rng) {
  f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/violin/Violin.",genes.sub[k],".pdf")
  print(paste0(k,": ",f.out))
  pdf(f.out,width=10,height=5)
  print(VlnPlot(dfseurat.sub,features=genes.sub[k],sort=FALSE) + NoLegend() + theme(axis.text.x=element_text(angle=90)))
  dev.off()
}




