library(Seurat)
library(Signac)

f = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds"
dfseurat = readRDS(f)
gene.activities <- GeneActivity(dfseurat)

dfseurat[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
dfseurat <- NormalizeData(
  object = dfseurat,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = 10000
)

f.out = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.GeneActivity.rds"
saveRDS(gene.activities,f.out)

library(data.table)
f.out = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.meta.txt"
fwrite(dfseurat@meta.data,f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)

f.out = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.GeneActivity.rds"
saveRDS(dfseurat,f.out)

df.sub = subset(dfseurat,subclust_v6=="HSCs")
Idents(df.sub) = "disease"
DefaultAssay(df.sub) = "GeneActivity"
(FindMarkers(df.sub, ident.1 = "T21", ident.2 = "H", test.use = "wilcox",features="GATA1",logfc.threshold=0,min.pct=0))
DefaultAssay(df.sub) = "chromvar"
(FindMarkers(df.sub, ident.1 = "T21", ident.2 = "H", test.use = "wilcox",slot='data',features="MA0035.4",logfc.threshold=0,min.pct=0,min.cells.feature=0,min.cells.group=0))

y1=as.numeric(subset(df.sub,disease=="T21")@assays$chromvar@data["MA0035.4",])
y2=as.numeric(subset(df.sub,disease=="H")@assays$chromvar@data["MA0035.4",])
wilcox.test(y1,y2)
t.test(y1,y2)
median(y1);median(y2)

df.sub@assays$chromvar["MA0035.4",]

print("getMatrixByID...")
library("JASPAR2022")
library("TFBSTools")
genelst <- getMatrixByID(JASPAR2022, ID = unique(rownames(df.sub@assays$chromvar)))

print("motif_name...")
motif_name <- as.character(
  unlist(
    mclapply(
      genelst,
      function(x) ID(x),
      mc.cores = 4)
  )
)

print("gene_name...")
gene_name <- as.character(
  unlist(
    mclapply(
      genelst,
      function(x) name(x),
      mc.cores = 4)
  )
)

print("motif-to-gene linking...")
res.all$i <- 1:nrow(res.all)
res.all.2 <- merge(res.all,
                   data.frame(motif_name,gene_name),
                   by.x='gene',
                   by.y='motif_name')

print("reordering...")
res.all.2 <- res.all.2[order(res.all.2$i),]




