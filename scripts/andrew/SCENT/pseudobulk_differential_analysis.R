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
# varfeat = FindVariableFeatures(df.aggre.norm)
# varfeat = varfeat[order(varfeat$vst.variance.standardized,decreasing = T),]
# varfeat.keep = rownames(varfeat)[1:2000]
# pca = prcomp(t(df.aggre.norm[varfeat.keep,]),scale.=TRUE)$x
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

cor.test(da_peaks.mg$avg_log2FC,da_peaks.mg$logFC)

df.sub2 = merge(df.sub,pb_da_peaks,by.x='peak',by.y='names',all.x=T)
df.sub2 = merge(df.sub2,de_genes,by.x='gene',by.y='names',all.x=T)
df.sub2 = merge(df.sub2,res.df,by.x='gene',by.y='names',all.x=T)
# subset(df.sub2,gene=="TFR2")
subset(df.sub2,peak=="chr14-65042220-65043484")
sum(subset(df.sub2,!duplicated(peak))$P.Value.x < 0.05)
sum(subset(df.sub2,!duplicated(gene))$P.Value.y < 0.05)
sum(subset(df.sub2,!duplicated(gene))$P.Value < 0.05)
sum(df.sub2$P.Value.x < 0.05 | df.sub2$P.Value.y < 0.05 | df.sub2$P.Value < 0.05)
sum(df.sub2$P.Value.x < 0.05 | df.sub2$P.Value.y < 0.05)
sum(df.sub2$P.Value.x < 0.05 | df.sub2$P.Value < 0.05)
sum(df.sub2$P.Value.x < 0.05 & df.sub2$P.Value.y < 0.05)
df.sub2
mean(df.sub2$)

subset(df.sub2,fdr.y < 0.01)

library(meta)

mean(subset(de_res.sub,de_res.sub$names %in% res.df.mg$gene)$P.Value.x < 0.05)
mean(subset(de_res.sub,de_res.sub$names %in% res.df.mg$gene)$P.Value.y < 0.05)
mean(subset(de_res.sub,de_res.sub$names %in% res.df.mg$gene)$P.Value.x < 0.05 | subset(de_res.sub,de_res.sub$names %in% res.df.mg$gene)$P.Value.y < 0.05)
mean(subset(pb_da_peaks,pb_da_peaks$names %in% res.df.mg$peak)$P.Value < 0.05)

mean(de_res.sub$P.Value.x < 0.05 | de_res.sub$P.Value.y < 0.05)
mean(pb_da_peaks)

sum(df.sub2$P.Value.x < 0.05)
sum(df.sub2$P.Value.y < 0.05)
sum(df.sub2$P.Value < 0.05)
sum(df.sub2$P.Value.y < 0.05 | df.sub2$P.Value < 0.05)

mean(df.sub2$P.Value.x < 0.05 | df.sub2$P.Value.y < 0.05)
df.sub2[!(df.sub2$P.Value.x < 0.05 | df.sub2$P.Value.y < 0.05),]
df.sub2[!(df.sub2$P.Value.x < 0.05 | df.sub2$P.Value < 0.05),]
df.sub2[!(df.sub2$P.Value.x < 0.05 | (df.sub2$P.Value < 0.05 | df.sub2$P.Value.y < 0.05)),]

df.sub2[(df.sub2$P.Value.x < 0.05),]


mean(df.sub2$P.Value.x > 0.05 & df.sub2$P.Value.y > 0.05)

subset(pb_da_peaks,names %in% df.sub$peak)
subset(res.df,names %in% df.sub$gene)


de_res = merge(res.df,de_genes,by.x='names',by.y='names',all=T)
cor.test(de_res$logFC.x,de_res$logFC.y)$estimate
cor.test(de_res$logFC.x,de_res$logFC.y)$p.value
ggplot(de_res,aes(x=logFC.x,y=logFC.y)) + geom_point() + geom_smooth(method='lm') + labs(x="logFC small multiome",y="logFC large scRNA") + theme_bw()
de_res.sub = subset(de_res,!is.na(logFC.x) & !is.na(logFC.y))
de_res.sub[order(de_res.sub$logFC.y,decreasing = T)[1:10],]
de_res.sub[order(de_res.sub$logFC.x,decreasing = F)[1:10],]

mean(de_res.sub$P.Value.x<0.05,na.rm=T)
mean(de_res.sub$P.Value.y<0.05,na.rm=T)


