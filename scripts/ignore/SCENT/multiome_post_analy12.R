dfseurat=readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

peak_of_interest="chr6-135096824-135097623"
gene_of_interest = "HBS1L"

dfseurat@assays$ATAC@counts[this_peak,]

DefaultAssay(dfseurat) <- "ATAC"
# dfseurat = NormalizeData(dfseurat,normalization.method = "LogNormalize", scale.factor = 10000)
dfseurat <- RunTFIDF(dfseurat)
Idents(dfseurat) = "disease"
da_peaks = FindMarkers(dfseurat,ident.1="T21",ident.2="H",test.use="wilcox",min.pct=0.05,logfc.threshold=0.05)
da_peaks$id = rownames(da_peaks)
# da_peaks = FindMarkers(dfseurat,ident.1="T21",ident.2="H",test.use="LR",min.pct=0.05,logfc.threshold=0.05)
fwrite(da_peaks,"/home/amarder/tmp/da_peaks.hsc.txt",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

DefaultAssay(dfseurat) <- "RNA"
library(Matrix.utils)
library(edgeR)
df.aggre <- aggregate.Matrix(
  t(
    dfseurat@assays$RNA@counts
    # dfseurat@assays$ATAC@counts
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
fwrite(pca.aggre,"/home/amarder/tmp/pca.hsc.txt",quote = F,na = "NA",row.names = F,col.names = T)
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
fwrite(res.df,"/home/amarder/tmp/pb_de.hsc.txt",row.names = F,col.names = T,quote = F,na = "NA",sep = "\t")

de_genes = fread("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
res.df.mg = merge(res.df,de_genes,by='names')
cor.test(res.df.mg$logFC.x,res.df.mg$logFC.y)
cor.test(res.df.mg$logFC.x,res.df.mg$logFC.y)$p.value
mean(res.df.mg$adj.P.Val.x < 0.1); mean(res.df.mg$adj.P.Val.y < 0.1)
res.df.mg[order(res.df.mg$adj.P.Val.y)[1:5],]
res.df.mg = subset(res.df.mg,)
df.sub2 = merge(df.sub2,de_genes,by.x='gene',by.y='names')

res.df= fread("~/Downloads/pb_de.hsc.txt",data.table = F,stringsAsFactors = F)
subset(res.df,names %in% df.sub$gene)
da_peaks = fread("~/Downloads/da_peaks.hsc.txt",data.table = F,stringsAsFactors = F)
de_genes = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
df.sub2 = merge(df.sub,da_peaks,by.x='peak',by.y='id',all.x=T)
df.sub2 = merge(df.sub2,de_genes,by.x='gene',by.y='names',all.x=T)
df.sub2 = merge(df.sub2,res.df,by.x='gene',by.y='names',all.x=T)
ggplot(df.sub2,aes(x=logFC.x,y=logFC.y)) + geom_point()

dfseurat.sub <- NormalizeData(dfseurat.sub,normalization.method = "LogNormalize", scale.factor = 10000)
  
trait_peaks = fread("/home/amarder/tmp/rbc_peaks_overlap.bed",data.table = F,stringsAsFactors = F)
snp_peaks = unique(subset(trait_peaks,V5 > 0.2)$V9)
length(snp_peaks)
nrow(da_peaks)
sum(da_peaks$p_val_adj < 0.1)
sum(subset(da_peaks,id %in% snp_peaks)$p_val_adj < 0.1)
snp_peaks = unique(subset(trait_peaks,V5 > 0.2)$V9)


da_peaks(subset(trait_peaks,V5 > 0.2))
FindMarkers(dfseurat,ident.1="T21",ident.2="H",test.use="wilcox",min.pct=0.05,logfc.threshold=0,
            features="chr6-135096824-135097623")

DefaultAssay(dfseurat) <- "RNA"
kclust=3;kmeans_results <- kmeans(Embeddings(dfseurat,reduction="harmony"),kclust)
data.frame()
dfseurat.sub.rna = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA.HSC_only.rds")
kmeans_results <- kmeans(Embeddings(dfseurat.sub.rna,reduction="rna_harmony"),kclust)
tmp1=data.frame(cell_dataset=dfseurat.sub.rna$cell_dataset,kmeans=as.factor(kmeans_results$cluster))
tmp<-data.frame(cell_dataset=dfseurat$cell_dataset,
  cell=colnames(dfseurat@assays$ATAC@counts),
  dataset=dfseurat$dataset,
                        disease=dfseurat$disease,
                        atac=as.numeric(dfseurat@assays$ATAC@counts[peak_of_interest,]>0),
                        atac_v2 = as.numeric(dfseurat@assays$ATAC@counts[peak_of_interest,]),
                        rna=as.numeric(dfseurat@assays$RNA@counts[gene_of_interest,]),
                        rna_v2=as.numeric(dfseurat@assays$RNA@data[gene_of_interest,])
                        )
tmp = merge(tmp,tmp1,by="cell_dataset")
aggregate(.~disease,data=tmp[,-c(1:3)],mean)
aggregate(.~disease+kmeans,data=tmp[,-c(1:3)],mean)
fisher.test(tmp$atac,tmp$disease)$p.value
# wilcox.test(tmp$atac_v2[tmp$disease=="H"],tmp$atac_v2[tmp$disease=="T21"])$p.value
table(tmp$kmeans,tmp$disease)
table(tmp$kmeans,tmp$dataset)
table(tmp$dataset,tmp$atac,tmp$disease)
table(tmp$kmeans,tmp$atac)

aggregate(.~disease+kmeans,data=tmp[,-c(1:3,5:7)],median)
aggregate(.~disease+kmeans+atac,data=tmp[,-c(1:3,6:7)],median)
summary(lm(rna_v2~disease+kmeans+atac,data=tmp[,-c(1:3,6:7)]))
fwrite(tmp,"~/tmp/HBS1L.expression.txt",quote=F,na="NA",row.names=F,col.names=T,sep = "\t")

tmp = fread("~/Downloads/HBS1L.expression.txt",data.table = F,stringsAsFactors = F)
tmp$kmeans = as.factor(tmp$kmeans)
tmp$atac = as.factor(tmp$atac)
tmp$lbl = paste("Cluster",tmp$kmeans,"-",ifelse(tmp$atac==1,"Accessible","Not accessible"))
ggplot(tmp,aes(x=disease,y=rna_v2,fill=atac)) + geom_boxplot(outlier.shape = NA) + 
  theme_bw() + labs(y='Expression') + scale_fill_manual(labels=c("Not accessible","Accessible"),values=c("darkgreen","lightblue"))
ggplot(tmp,aes(x=disease,y=rna_v2,fill=lbl)) + geom_boxplot(outlier.shape = NA) + 
  theme_bw() + labs(y='Expression')

ggplot(tmp,aes(x=disease,y=rna_v2,fill=atac)) + geom_violin() + theme_bw()


wilcox.test(tmp$rna_v2[tmp$disease=="H"],tmp$rna_v2[tmp$disease=="T21"])$p.value
wilcox.test(tmp$rna_v2[tmp$disease=="T21" & tmp$kmeans==2],tmp$rna_v2[tmp$disease=="T21" & tmp$kmeans==1])$p.value

t.test(tmp$rna_v2[tmp$disease=="H"],tmp$rna_v2[tmp$disease=="T21"])$p.value

mean(tmp$rna[tmp$disease=="H"])
mean(tmp$rna[tmp$disease=="T21"])