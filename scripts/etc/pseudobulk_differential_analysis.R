DefaultAssay(dfseurat) <- "ATAC"
library(Matrix.utils)
library(edgeR)
df.aggre <- aggregate.Matrix(
  t(
    # dfseurat@assays$RNA@counts
    dfseurat@assays$ATAC@counts
  ),
  groupings=dfseurat$dataset,fun="sum")
df.aggre <- t(df.aggre)
df.aggre <- as.data.frame(df.aggre)
df.aggre.norm = cpm(df.aggre)
pca = prcomp(t(df.aggre.norm[rowSums(df.aggre.norm)>0,]),scale.=TRUE)$x
pca.aggre = merge(
  unique(dfseurat@meta.data[,c("dataset","disease")]),
  pca[,1:2],
  by.x="dataset",
  by.y="row.names"
)
fwrite(pca.aggre,"/home/amarder/tmp/pca_atac.hsc.txt",quote = F,na = "NA",row.names = F,col.names = T)
pca.aggre = fread("~/Downloads/pca.hsc.txt",data.table = F,stringsAsFactors = F)
ggplot(pca.aggre,aes(x=PC1,y=PC2,col=disease)) + geom_point() + theme_bw() + labs(title="HSC Pseudobulks") + theme(plot.title = element_text(hjust=0.5))

metadata_to_use = unique(dfseurat@meta.data[,c("dataset","disease")])
rownames(metadata_to_use) = NULL
metadata_to_use = metadata_to_use[match(colnames(df.aggre),metadata_to_use$dataset),]
rownames(metadata_to_use) = metadata_to_use$dataset

# Standard usage of limma/voom
geneExpr = DGEList( df.aggre )
keep <- filterByExpr(geneExpr, group=metadata_to_use$disease)
geneExpr <- geneExpr[keep,]

geneExpr = calcNormFactors( geneExpr )

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel
param = SnowParam(8, "SOCK", progressbar=TRUE)

# The variable to be tested must be a fixed effect
form <- ~ disease  

# A positive FC is increased expression in the DS compared to healthy
metadata_to_use$disease <- factor(metadata_to_use$disease,levels=c("H","T21"))

# estimate weights using linear mixed model of dream
library(variancePartition); library('BiocParallel')
vobjDream = voomWithDreamWeights( geneExpr, form, metadata_to_use, BPPARAM=param )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream( vobjDream, form, metadata_to_use )
fitmm = eBayes(fitmm)

# reorganize:
# res.df <- as.data.frame(topTable( fitmm, number=Inf))
res.df <- as.data.frame(topTable( fitmm, number=Inf,coef = "diseaseT21" ))
res.df$names <- rownames(res.df)
fwrite(res.df,"/home/amarder/tmp/pb_de_atac.hsc.txt",row.names = F,col.names = T,quote = F,na = "NA",sep = "\t")

pca.aggre = fread("~/Downloads/pca_atac.hsc.txt",data.table = F,stringsAsFactors = F)
ggplot(pca.aggre,aes(x=PC1,y=PC2,col=disease)) + geom_point() + theme_bw() + labs(title="HSC Pseudobulks") + theme(plot.title = element_text(hjust=0.5))

da_peaks = fread("~/Downloads/da_peaks.hsc.txt",data.table = F,stringsAsFactors = F)
pb_da_peaks = fread("~/Downloads/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)
da_peaks.mg = merge(da_peaks,pb_da_peaks,by.x="id",by.y="names")
ggplot(da_peaks.mg,aes(x=avg_log2FC,y=logFC)) + geom_point() + geom_abline(slope=1,intercept=0,col='red',lty='dashed') + theme_bw()
ggplot(da_peaks.mg,aes(x=-log10(p_val),y=-log10(P.Value))) + geom_point() + geom_abline(slope=1,intercept=0,col='red',lty='dashed') + theme_bw()


