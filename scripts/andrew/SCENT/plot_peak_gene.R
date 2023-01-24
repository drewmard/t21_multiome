library(data.table)

# 0. SCENT!
celltype_to_use="HSCs_H"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out2/",celltype_to_use,".txt")
res.df.h = fread(output_file,data.table = F,stringsAsFactors = F)

celltype_to_use="HSCs_T21"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out2/",celltype_to_use,".txt")
res.df.t21 = fread(output_file,data.table = F,stringsAsFactors = F)

# ggplot(res.df.t21,aes(x=-log10(p),y=-log10(pval))) + geom_point() + geom_abline(slope=1,intercept=0)

res.df.mg = merge(res.df.h,res.df.t21,all=TRUE,by=c("gene","peak"))

cor(res.df.mg$beta.x,res.df.mg$beta.y,use='na.or.complete')
tmp = subset(res.df.mg,fdr.x < 0.2 & fdr.y < 0.2)
table(sign(tmp$beta.x)!=sign(tmp$beta.y))

subset(res.df.mg,fdr.x < 0.2 & fdr.y < 0.2 & sign(beta.x)!=sign(beta.y))

##########################################################################################
##########################################################################################

library(data.table)
library(ggplot2)
library(Seurat)
library(Signac)
dfseurat = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

gene_of_interest = "RAPGEF2"; peak_of_interest = "chr4-159150290-159152052"

gene_ct = dfseurat@assays$RNA@counts[gene_of_interest,]
gene=dfseurat@assays$RNA@data[gene_of_interest,]
this_peak <- dfseurat@assays$ATAC@data[peak_of_interest,]
this_peak.bin <- dfseurat@assays$ATAC@counts[peak_of_interest,]
this_peak.bin = this_peak.bin
this_peak.bin[this_peak.bin > 1] <- 1
df <- data.frame(rna=gene,rna_raw=gene_ct,atac=this_peak,atac.bin = this_peak.bin,disease=dfseurat@meta.data$disease)
f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/plots/peak_gene_link.",gene_of_interest,".",peak_of_interest,".txt")
fwrite(df,file=f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
aggregate(mrna~atac,df,median)
aggregate(mrna~atac,subset(df,disease=="H"),median)
aggregate(mrna~atac,subset(df,disease=="T21"),median)
aggregate(mrna~atac,subset(df,disease=="H"),length)
aggregate(mrna~atac,subset(df,disease=="T21"),length)

#########

library(data.table)
library(ggplot2)
gene_of_interest = "RAPGEF2"; peak_of_interest = "chr4-159150290-159152052"
f = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/plots/peak_gene_link.",gene_of_interest,".",peak_of_interest,".txt")
df = fread(f,data.table = F,stringsAsFactors = F)
head(df)
ggplot(df,aes(x=disease,y=rna,fill=as.factor(atac.bin))) + geom_boxplot() + 
  labs(fill="chr4-\n159150290-\n159152052\n\nAccessibility",x='Disease status',y="RAPGEF2 Expression") +
  scale_fill_manual(values=c('lightblue','orange'),labels=c("Not open","Open")) +
  theme_bw() + theme(panel.grid = element_blank())
# ggplot(df,aes(x=atac,y=rna,col=as.factor(disease))) + geom_point()






