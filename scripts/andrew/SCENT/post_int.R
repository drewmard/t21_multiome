library(data.table)
df = fread("~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs.DE.txt",data.table = F,stringsAsFactors = F)
# subset(df,gene=="AFF3" & peak=="chr2-100141820-100142865")

library(data.table)
df = fread("~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs.DE.txt",data.table = F,stringsAsFactors = F)

dfint = fread("~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs_all.all.txt",data.table = F,stringsAsFactors = F)
df.mg = merge(df,dfint,by = c("gene",'peak'))

df$gene_peak = paste(df$peak,df$gene)
dfint$gene_peak = paste(dfint$peak,dfint$gene)
sum(!(df$gene_peak %in% dfint$gene_peak))
sum(duplicated(df$gene_peak))
df[duplicated(df$gene_peak),][1:3,]
subset(df,gene=="HERC3" & peak=="chr4-88612949-88613880")
sum(duplicated(df$gene_peak))

df[7973,]

df.mg.sub = subset(df.mg,fdr_H < 0.2 | fdr_t21 < 0.2)
df.mg.sub$fdr = p.adjust(df.mg.sub$pval,method='fdr')

sum(df.mg$fdr_H < 0.2 & df.mg$fdr_t21 < 0.2,na.rm = T)

mean(df.mg$fdr < 0.2)
plot(df.mg$beta_t21-df.mg$beta_H,df.mg$beta)
nrow(df.mg.sub)
sum(df.mg.sub$fdr < 0.2)
mean(df.mg.sub$fdr < 0.2)

mean(subset(df.mg.sub,fdr_H < 0.2 & fdr_t21 < 0.2)$fdr < 0.2)
mean(subset(df.mg.sub,fdr_H < 0.2 & (fdr_t21 > 0.2 | is.na(fdr_t21)))$fdr < 0.2)
mean(subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H)))$fdr < 0.2)
nrow(subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H))))

H_only = unique(subset(df.mg.sub,fdr_H < 0.2 & (fdr_t21 > 0.2 | is.na(fdr_t21)) & fdr < 0.2)$gene)
t21_only = unique(subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H)) & fdr < 0.2)$gene)
t21_only_but_still_access_and_expressed = unique(subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2) & fdr < 0.2)$gene)
disconcordant = unique(subset(df.mg.sub,fdr_H < 0.2 & fdr_t21 < 0.2 & fdr < 0.2 & sign(beta_H)!=sign(beta_t21))$gene)
t21_enhanced_but_still_access_and_expressed = unique(subset(df.mg.sub,fdr_H < 0.2 & fdr_t21 < 0.2 & fdr < 0.2 & sign(beta_H)==sign(beta_t21) & sign(beta)==sign(beta_t21))$gene)
t21_diminished_but_still_access_and_expressed = unique(subset(df.mg.sub,fdr_H < 0.2 & fdr_t21 < 0.2 & fdr < 0.2 & sign(beta_H)==sign(beta_t21) & sign(beta)!=sign(beta_t21))$gene)
t21_enhanced = unique(subset(df.mg.sub,(fdr_H < 0.2 | is.na(fdr_H)) & fdr_t21 < 0.2 & fdr < 0.2 & !(gene %in% disconcordant) & sign(beta)==sign(beta_t21))$gene)
t21_diminished = unique(subset(df.mg.sub,(fdr_H < 0.2 | is.na(fdr_H)) & fdr_t21 < 0.2 & fdr < 0.2 & !(gene %in% disconcordant) & sign(beta)!=sign(beta_t21))$gene)
any_interaction = unique(subset(df.mg.sub,fdr < 0.2)$gene)

fwrite(data.frame(H_only),"~/Downloads/H_only",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(t21_only),"~/Downloads/t21_only",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(t21_only_but_still_access_and_expressed),"~/Downloads/t21_only_but_still_access_and_expressed",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(disconcordant),"~/Downloads/disconcordant",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(t21_enhanced_but_still_access_and_expressed),"~/Downloads/t21_enhanced_but_still_access_and_expressed",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(t21_diminished_but_still_access_and_expressed),"~/Downloads/t21_diminished_but_still_access_and_expressed",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(t21_enhanced),"~/Downloads/t21_enhanced",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(t21_diminished),"~/Downloads/t21_diminished",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
fwrite(data.frame(any_interaction),"~/Downloads/any_interaction",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

t.test(subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H)) & fdr > 0.2)$logFC_lgRNA,na.rm=T)
# tmp = subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H)) & fdr > 0.2)
tmp = subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H)) & fdr > 0.2 & !duplicated(gene))
# tmp$gene
# mean(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0,na.rm=T)
tab = table(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0)
large_RNA_pb = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
mu = mean(large_RNA_pb$P.Value < 0.05 & large_RNA_pb$logFC > 0,na.rm=T)
binom.test(tab["TRUE"],tab["TRUE"] + tab["FALSE"],mu)

tab = table(tmp$adj.P.Val_lgRNA < 0.05 & tmp$logFC_lgRNA > 0)
mu = mean(large_RNA_pb$adj.P.Val < 0.05 & large_RNA_pb$logFC > 0,na.rm=T)
binom.test(tab["TRUE"],tab["TRUE"] + tab["FALSE"],mu)

# t21 specific interactions
tmp = subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H)) & fdr < 0.2)
mean(sign(tmp$beta)==sign(tmp$beta_t21))
tmp = subset(df.mg.sub,fdr_H < 0.2 & (fdr_t21 > 0.2 | is.na(fdr_t21)) & fdr < 0.2)
mean(sign(tmp$beta)!=sign(tmp$beta_H))
tmp = subset(df.mg.sub,fdr_t21 < 0.2 & fdr_H < 0.2 & fdr < 0.2 & gene != "RAPGEF2")
mean(sign(tmp$beta)==sign(tmp$beta_t21))
tmp$enhanced = as.numeric(sign(tmp$beta)==sign(tmp$beta_t21))
mean(tmp$enhanced)

tmp = subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2) & fdr < 0.2)
tmp$enhanced = as.numeric(sign(tmp$beta)==sign(tmp$beta_t21))
tmp[,c("beta","beta_H","beta_t21","enhanced","P.Value_lgRNA")]
median(tmp$P.Value_lgRNA,na.rm = T)
median(tmp$logFC_lgRNA,na.rm = T)

tmp["26312",]

tmp = subset(df.mg.sub,fdr_t21 < 0.2 & (fdr_H > 0.2 | is.na(fdr_H)) & fdr < 0.2)
mean(sign(tmp$beta)==sign(tmp$beta_t21))


tmp = subset(df.mg.sub,fdr_H < 0.2 & (fdr_t21 > 0.2 | is.na(fdr_t21)) & fdr < 0.2)
mean(sign(tmp$beta)!=sign(tmp$beta_H))

subset(df.mg,fdr_H < 0.2 & fdr_t21 < 0.2 & sign(beta_H)!=sign(beta_t21))

subset(df,gene=="CRIM1" & peak=="chr2-36355063-36358149")
subset(df,gene=="TFR2" & peak=="chr7-100639909-100642992")
