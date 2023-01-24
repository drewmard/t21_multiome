library(Seurat)
library(Signac)
library(Matrix.utils)
library(edgeR)
library(DESeq2)
library('variancePartition')
library('BiocParallel')

# Data input:
dfseurat=readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds")
# dfseurat=readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

annot <- getBM(attributes = c('hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol',
               values = rownames(df.aggre),
               mart = ensembl)
annot2 = subset(annot,chromosome_name %in% seq(1,22))
rownames(dfaggre) %in% annot2$hgnc_symbol

df1 <- merge(res.df1,annot[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
df1 <- df1[df1$chromosome_name %in% seq(1,22),]
df1$chromosome_name <- factor(df1$chromosome_name,levels=seq(1,22))
df1$chr21 <- factor(ifelse(df1$chromosome_name==21,'Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))


dataset.uniq = unique(dfseurat$dataset)
# celltype.lst=c("HSCs","MEMPs","Early erythroid","Late erythroid")
celltype.lst=c("HSCs","MEMPs","Early erythroid")

tab = table(dfseurat$dataset,dfseurat$subclust_v6)
cellSize = apply(tab[,celltype.lst],1,min)

i=1
j=1

for (i in 1:length(celltype.lst)) { 
  celltype=celltype.lst[i]
  ind.lst=list()
  for (j in 1:length(dataset.uniq)) { 
    datasetUse = dataset.uniq[j]
    which.cell = dfseurat$subclust_v6==celltype
    which.data = dfseurat$dataset==datasetUse
    set.seed(031995)
    # set.seed(03191995)
    ind.all_cell = which(which.cell & which.data)
    ind.lst[[datasetUse]] = sample(ind.all_cell,cellSize[datasetUse])
  }
  ind = unlist(ind.lst)
  dfsub = dfseurat[,ind]
  
  
  # dfsub = dfseurat[,which(dfseurat$subclust_v6=="HSCs")]
  # Idents(dfsub) = 'disease'
  # x = FindMarkers(dfsub,assay="RNA",ident.1="T21",ident.2="H")
  # df1 <- merge(x,annot2[,c("hgnc_symbol","chromosome_name","start_position")],by.x="row.names",by.y="hgnc_symbol")
  # aggregate(df1$avg_log2FC,by=list(df1$chromosome_name==21),mean)
  
  df.aggre <- aggregate.Matrix(
    t(
      # dfseurat@RNA$counts[,ind]
      GetAssayData(object = dfsub, slot = "counts", assay="RNA")
    ),
    groupings=dfsub$dataset,fun="sum")
  
  df.aggre <- t(df.aggre)
  df.aggre <- as.data.frame(df.aggre)
  
  x <- unique(dfsub@meta.data[,c("dataset","disease")])
  rownames(x) <- x[,"dataset"]
  
  metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
  df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
  
  # Standard usage of limma/voom
  geneExpr = DGEList( df.aggre )
  keep <- filterByExpr(geneExpr, group=metadata_to_use$disease)
  geneExpr <- geneExpr[keep,]
  # geneExpr = DGEList( df.aggre[rownames(df.aggre) %in% annot2$hgnc_symbol,] )
  geneExpr = calcNormFactors( geneExpr )
  
  ###########
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(8, "SOCK", progressbar=TRUE)
  
  # The variable to be tested must be a fixed effect
  form <- ~ disease
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )
  
  # Fit the model on each gene
  fitmm = dream( vobjDream, form, metadata_to_use )
  fitmm = eBayes(fitmm)
  
  # reorganize:
  res.df1 <- as.data.frame(topTable( fitmm, number=Inf,coef = "diseaseT21" ))
  res.df1$names <- rownames(res.df1)
  
  df1 <- merge(res.df1,annot2[,c("hgnc_symbol","chromosome_name","start_position")],by.x="names",by.y="hgnc_symbol")
  # df1$chr21 <- factor(ifelse(df1$chromosome_name=='21','Chr 21','Not Chr 21'),levels=c('Not Chr 21','Chr 21'))
  
  print(aggregate(df1$logFC,by=list(df1$chromosome_name==21),mean))
  print(aggregate(df1$adj.P.Val < 0.05,by=list(chr21=df1$chromosome_name==21,upreg=df1$logFC > 0),sum))
  
  library(data.table)
  fwrite(df1,paste0("~/tmp/",celltype,".txt"),quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
}
# aggregate(df1$AveExpr,by=list(chr21=df1$chromosome_name==21,upreg=df1$logFC > 0),mean)



library(data.table)
for (i in 1:length(celltype.lst)) { 
  celltype=celltype.lst[i]
df1 = fread(paste0("~/tmp/",celltype,".txt"),data.table = F,stringsAsFactors = F)
print(celltype)
print("Avg LFC:")
print(aggregate(df1$logFC,by=list(chr21=df1$chromosome_name==21),mean))
print("Number of DEGs:")
print(aggregate(df1$adj.P.Val < 0.05,by=list(chr21=df1$chromosome_name==21,upreg=df1$logFC > 0),sum))
}


library(data.table)
for (i in 1:length(celltype.lst)) { 
  celltype=celltype.lst[i]
  df1 = fread(paste0("~/tmp/",celltype,".txt"),data.table = F,stringsAsFactors = F)
  print(celltype)
  # print(subset(df1,names=="LINC01478")[,1:6])
  # print(sum(df1[grep("LINC",df1$names),6] < 0.05))
  print(length(df1[grep("LINC",df1$names),"logFC"]))
  # print(mean(df1[grep("LINC",df1$names),"logFC"]))
}

library(data.table)
df1 = list()
for (i in 1:length(celltype.lst)) { 
  celltype=celltype.lst[i]
  df1[[celltype]] = fread(paste0("~/tmp/",celltype,".txt"),data.table = F,stringsAsFactors = F)
  print(celltype)
  if (i==1) {
    dfall = df1[[celltype]]
  } else {
    dfall = merge(dfall,df1[[celltype]],by='names',all=T)
  }
}

wilcox.test(dfall$logFC.x[grep("LINC",dfall$names)],
            dfall$logFC[grep("LINC",dfall$names)])
dfall[grep("LINC",dfall$names),][1:2,]

tmp = dfall[grep("LINC",dfall$names),]
tmp[order(tmp$adj.P.Val)[1:3],]
aggregate(tmp$adj.P.Val.x < 0.05,by=list(upreg=tmp$logFC.x > 0),sum)
aggregate(tmp$adj.P.Val.y < 0.05,by=list(upreg=tmp$logFC.y > 0),sum)


