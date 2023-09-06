# see: cycHSC_vs_lesscycHSC_multiome
DefaultAssay(dfseurat) <- "ATAC"
library(Matrix.utils)
library(edgeR)

k=2
dfseurat$clust = as.numeric(dfseurat$kmeans_RNA==k) + 1

df.aggre.lst=list()
metadata_to_use.lst =list()
# for (k in 1:2) {
for (k in 1:2) {
  rm(df.sub)
  df.sub=dfseurat[,dfseurat$clust==k & dfseurat$disease=="T21"]
  # df.sub=dfseurat[,dfseurat$lin==k & dfseurat$disease=="T21"]
  
  df.aggre <- aggregate.Matrix(
    t(
      df.sub@assays$RNA@counts
      # df.sub@assays$ATAC@counts
    ),
    groupings=df.sub$dataset,fun="sum")
  df.aggre <- t(df.aggre)
  df.aggre <- as.data.frame(df.aggre)
  
  metadata_to_use = unique(df.sub@meta.data[,c("dataset","clust")])
  # metadata_to_use = unique(df.sub@meta.data[,c("dataset","lin")])
  
  df.aggre.lst[[k]] = df.aggre
  metadata_to_use.lst[[k]] = metadata_to_use
  
  # for (k in 1:2) { 
  colnames(df.aggre.lst[[k]]) = paste0(colnames(df.aggre.lst[[k]])," ",k)
  metadata_to_use.lst[[k]]$dataset = paste0(metadata_to_use.lst[[k]]$dataset," ",k)
  # }
}

# df.aggre = do.call(cbind,df.aggre.lst[2:3])
# metadata_to_use = do.call(rbind,metadata_to_use.lst[2:3])

df.aggre = do.call(cbind,df.aggre.lst)
metadata_to_use = do.call(rbind,metadata_to_use.lst)

metadata_to_use=metadata_to_use[match(colnames(df.aggre),metadata_to_use$dataset),]
df.aggre.norm = cpm(df.aggre)
df.aggre.norm["MKI67",]

# Standard usage of limma/voom
geneExpr = DGEList( df.aggre )
# keep <- filterByExpr(geneExpr, group=metadata_to_use$kmeans_RNA)
keep <- filterByExpr(geneExpr, group=metadata_to_use$clust)
geneExpr <- geneExpr[keep,]

geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ clust  

# A positive FC is increased expression in the DS compared to healthy
# metadata_to_use$kmeans_RNA <- factor(metadata_to_use$kmeans_RNA,levels=c("3","2"))
# metadata_to_use$kmeans_RNA <- factor(metadata_to_use$kmeans_RNA,levels=c("1","2"))
# metadata_to_use$lin <- factor(metadata_to_use$lin,levels=c("1","2"))
metadata_to_use$clust <- factor(metadata_to_use$clust,levels=c("1","2"))
# metadata_to_use$lin <- factor(metadata_to_use$lin,levels=c("3","2"))

# estimate weights using linear mixed model of dream
library(variancePartition); library('BiocParallel')
vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata_to_use )
fitmm = eBayes(fitmm)

# reorganize:
# res.df <- as.data.frame(topTable( fitmm, number=Inf))
res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "clust2" ))
# res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "kmeans_RNA2" ))
res.df$names <- rownames(res.df)

subset(res.df,names=="MKI67")


# fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

# scrna = fread("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/differential_analysis/HSC_v_CyclingHSC.txt",data.table = F,stringsAsFactors = F)
res.df.mg = merge(scrna,res.df,by='names')
cor.test(res.df.mg$logFC.x,res.df.mg$logFC.y)
res.df.mg.sub = subset(res.df.mg,adj.P.Val.y < 0.1)
cor.test(res.df.mg.sub$logFC.x,res.df.mg.sub$logFC.y)
fisher.test(res.df.mg.sub$logFC.x > 0,res.df.mg.sub$logFC.y > 0)



