library(Seurat)
library(Signac)
f="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds"
df = readRDS(f)
DefaultAssay(df) <- "RNA"

table(df$dataset)
genes_to_use = c("CD34","SPINK2","PROM1","MLLT3","TESPA1",
                 "GATA2","GATA1","ALAS2","HBA1","HDC",
                 "ITGA2B","MPO","AZU1","SPI1","FCN1","VCAN",
                 "CTSB","NKG7","PRF1","IL2RB",
                 "IL7R","PAX5","IGHM","IRF8","CLEC4C",
                 "CLEC10A","MKI67")
cellOrder = c("HSCs",
              "MEMPs",
              "Early erythroid",
              "Late erythroid",
              "Cycling erythroid",
              "Mast cells",
              "Megakaryocytes",
              "Granulocyte progenitors",
              "Neutrophils",
              "Pro-inflammatory macrophages",
              "Kupffer cells",
              "NK cells",
              "T cells",
              "B cells",
              "pDCs",
              "cDCs",
              "Stroma")#,
              # "Unknown",
              # "No marker")
# cellOrder <- rev(cellOrder)
df@meta.data$cluster_label <- factor(as.character(df@meta.data$subclust_v6),levels=cellOrder)

unique(subset(df@meta.data,is.na(cluster_label))$subclust_v6)
# unique(subset(df.sub@meta.data,is.na(cluster_label))$subclust_v6)

df.sub = subset(df,dataset=="15582 nuclei A")
Idents(df.sub) = "cluster_label"
df.sub = ScaleData(df.sub)
pdf("/oak/stanford/groups/smontgom/amarder/tmp/heatmap.15582_nuclei_A.pdf",width = 20,height=10)
print(DoHeatmap(df.sub,features = genes_to_use,raster=FALSE))
dev.off()
df.sub = subset(df,dataset=="15582 nuclei B")
Idents(df.sub) = "cluster_label"
df.sub = ScaleData(df.sub)
pdf("/oak/stanford/groups/smontgom/amarder/tmp/heatmap.15582_nuclei_B.pdf",width = 20,height=10)
print(DoHeatmap(df.sub,features = genes_to_use,raster = FALSE))
dev.off()

pdf("/oak/stanford/groups/smontgom/amarder/tmp/heatmap.15582_nuclei_B.sub.pdf",width = 20,height=10)
print(DoHeatmap(df.sub[,1:40],features = genes_to_use))
dev.off()

# x=40; png("/oak/stanford/groups/smontgom/amarder/tmp/heatmap.15582_nuclei_B.sub.png",width = 20*x,height=10*x)
# print(DoHeatmap(df.sub[,1:40],features = genes_to_use,raster = FALSE))
# dev.off()

df = fread("~/Documents/Research/t21-proj/out/full/data/meta.10X_Healthy_Liver.umap2d.cells_removed.txt",data.table = F,stringsAsFactors = F)
df.sub = subset(df,sample=="15582 B")


# mamba activate r
library(Seurat)

f = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_DownSyndrome_Liver.umap2d.cells_removed.rds"
df = readRDS(f)
df@assays$RNA@key <- "rna_"
# df.sub = df[,df@meta.data$sample=="T21 15582 B"]
# sum(df@meta.data$sample=="T21 15582 B")
df.sub = subset(df,sample=="T21 15582 B")
f.out = "/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/10X_DownSyndrome_Liver.umap2d.cells_removed.15582_B.rds"
saveRDS(df.sub,f.out)

# which(df$sample=="15582 T21 15582 B")

"leiden_v10"

df.sub = subset(df.sub,sample=="T21 15582 B" & leiden_v10 != "X")

genes_to_use = c("CD34","SPINK2","PROM1","MLLT3","TESPA1",
                 "GATA2","GATA1","ALAS2","HBA1","HDC",
                 "ITGA2B","MPO","AZU1","SPI1","FCN1","VCAN",
                 "CTSB","NKG7","PRF1","IL2RB",
                 "IL7R","PAX5","IGHM","IRF8","CLEC4C",
                 "CLEC10A","MKI67")
cellOrder = c("HSCs/MPPs",
              "Cycling HSCs/MPPs",
              "MEMPs",
              "Early erythroid cells",
              "Late erythroid cells",
              "Cycling erythroid",
              "Mast cells",
              "Megakaryocytes",
              "Cycling megakaryocytes",
              "Granulocyte progenitors",
              "Monocyte progenitors",
              "Neutrophils",
              "Inflammatory macrophages",
              "Kupffer cells",
              "NK progenitors",
              "NK cells",
              "T cells",
              "Pre pro B cells",
              "Pro B cells",
              "pDCs",
              "Cycling pDCs",
              "cDC2",
              "Stroma",
              "X")#,
# "Unknown",
# "No marker")
# cellOrder <- rev(cellOrder)
df.sub@meta.data$cluster_label <- factor(as.character(df.sub@meta.data$leiden_v10),levels=cellOrder)

unique(subset(df.sub@meta.data,is.na(cluster_label))$leiden_v10)
# unique(subset(df.sub@meta.data,is.na(cluster_label))$subclust_v6)

Idents(df.sub) = "cluster_label"
df.sub = ScaleData(df.sub)
pdf("/oak/stanford/groups/smontgom/amarder/tmp/heatmap.scrna.15582_B.pdf",width = 20,height=10)
print(DoHeatmap(df.sub,features = genes_to_use,raster=FALSE))
dev.off()



