
###################

library(data.table)
traitName="wbc"; res = list()
for(traitName in c("rbc","wbc","lymph")) {

f = paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt")
df = fread(f,data.table = F,stringsAsFactors = F)
df$trs_ex = df$TRS > quantile(df$TRS,probs=0.9)
df$trs_ex = ifelse(df$trs_ex,"GWAS-enriched HSCs","Other HSCs")

lineage = fread("~/Downloads/Lineage_Probabilities_per_Cell.csv",data.table = F,stringsAsFactors = F)
mapping = fread("~/Downloads/dataset_metadata_to_cell.csv",data.table = F,stringsAsFactors = F,header=T)
which(lineage$V1!=mapping$V1)
lineage = cbind(mapping,lineage[,-1])
lineage$cell = unlist(lapply(strsplit(lineage$V1,"_"),function(x) x[1]))
df2 = subset(df,subclust_v6=="HSCs")

df.mg = merge(df2,lineage,by=c("cell","dataset"),all=T)
df.mg$lin = apply(df.mg[,paste0("Lineage_",c(1:3),"_prob")],1,which.max)
aggregate(TRS~lin,df.mg,median)

tmp = aggregate(TRS~lin+dataset,df.mg,sum)
tmp2 = aggregate(TRS~lin+dataset,df.mg,length)
tmp = merge(tmp,tmp2,by=1:2)
tmp$lin = as.factor(tmp$lin)
contrasts(tmp$lin)
# Fit the ANCOVA model
tmp$lin1=tmp$lin=="1"
tmp$lin2=tmp$lin=="2"
tmp$lin3=tmp$lin=="3"
mod <- lm(TRS.x ~ lin1 + TRS.y, data = tmp)
x=summary(mod)$coef[2,4]
mod <- lm(TRS.x ~ lin2 + TRS.y, data = tmp)
y=summary(mod)$coef[2,4]
mod <- lm(TRS.x ~ lin3 + TRS.y, data = tmp)
z=summary(mod)$coef[2,4]
res[[traitName]] = data.frame(x,y,z)
}
res.df=melt(do.call(rbind,res))
res.df$fdr = p.adjust(res.df$value,method='fdr')
res.df

