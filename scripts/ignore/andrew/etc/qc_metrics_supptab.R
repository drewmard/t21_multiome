library(data.table)
type="ds"
f=paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/meta.",type,".txt")
df=fread(f,data.table = F,stringsAsFactors = F)
x1=aggregate(df[,c("nCount_ATAC","TSS.enrichment","nCount_RNA","nFeature_RNA")],by=list(df$dataset),median)
x2=aggregate(df[,c("nCount_ATAC")],by=list(df$dataset),length)
x = merge(x1,x2,by="Group.1")
x[c(5,6,1,2,3,4),]
fwrite(x[c(5,6,1,2,3,4),],"~/Downloads/tmp1.csv")

type="h"
f=paste0("/Users/andrewmarderstein/Documents/Research/neuro-variants/output/data/DS_Multiome_",type,"_v2/meta.",type,".txt")
df=fread(f,data.table = F,stringsAsFactors = F)
x1=aggregate(df[,c("nCount_ATAC","TSS.enrichment","nCount_RNA","nFeature_RNA")],by=list(df$dataset),median)
x2=aggregate(df[,c("nCount_ATAC")],by=list(df$dataset),length)
x = merge(x1,x2,by="Group.1")
# x
x[c("15828 C","15828 D"),]
fwrite(x,"~/Downloads/tmp1.csv")

