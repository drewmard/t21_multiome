library(data.table)

# 0. SCENT!
celltype_to_use="HSCs_H"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
res.df.h = fread(output_file,data.table = F,stringsAsFactors = F)

celltype_to_use="HSCs_T21"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out2/",celltype_to_use,".txt")
res.df.t21 = fread(output_file,data.table = F,stringsAsFactors = F)

# ggplot(res.df.t21,aes(x=-log10(p),y=-log10(pval))) + geom_point() + geom_abline(slope=1,intercept=0)

res.df.mg = merge(res.df.h,res.df.t21,all=TRUE,by=c("gene","peak"))

cor(res.df.mg$beta.x,res.df.mg$beta.y,use='na.or.complete')
subset(res.df.mg,fdr.x < 0.2 & fdr.y < 0.2 & sign(beta.x)!=sign(beta.y))



gene=chunkinfo$gene[i]
this_peak=chunkinfo$peak[i]
atac_target<-data.frame(cell=colnames(atac.all),atac=as.numeric(atac.all[this_peak,]))
mrna_target<-mrna[gene,]
df <- data.frame(cell=names(mrna_target),exprs=as.numeric(mrna_target))
df<-merge(df,atac_target,by="cell")
df<-merge(df,meta,by="cell")

# Subset cells to test:
df2 <- df[df$celltype==celltype_to_use,]

# Binarize peaks:
df2[df2$atac>0,]$atac<-1



tmp = subset(res.df.mg,fdr.x < 0.1 & fdr.y < 0.1); table(sign(tmp$beta.x) != sign(tmp$beta.y))
table(res.df.h$fdr < 0.2)
table(res.df.t21$fdr < 0.2)
table(res.df.mg$fdr.x < 0.2, res.df.mg$fdr.y < 0.2)
table(res.df.mg$fdr.x < 0.2 | res.df.mg$fdr.y < 0.2)
fisher.test(res.df.mg$fdr.x < 0.2, res.df.mg$fdr.y < 0.2)

# res.df.mg[order(res.df.mg$pval-res.df.mg$p.y,decreasing = T),][1,]
ggplot(res.df.mg,aes(x=beta.x,y=beta.y)) + geom_point() + geom_smooth(method='lm')
ggplot(res.df.mg,aes(x=-log10(fdr.x),y=-log10(fdr.y),col=sign(beta.x)!=sign(beta.y))) + geom_point() #+ geom_smooth(method='lm')


# res.df.t21[order(res.df.t21$p - res.df.t21$pval,decreasing = T),][1:5,]
# seq(0,1,length.out=nrow(res.df.t21))
# tmp = data.frame(exp=-log10(seq(0,1,length.out=nrow(res.df.t21))),
# obs=sort(-log10(res.df.t21$pval),decreasing = T))
# ggplot(tmp,aes(x=exp,y=obs)) +  theme_bw() +
#   geom_point() + geom_abline(slope=1,intercept=0)
# ggplot(tmp,aes(x=exp,y=obs)) +  theme_bw() +
#   geom_point() + geom_abline(slope=1,intercept=0)

traitName="rbc"
# traitName="lymph"


# 1. find trait-associated SNPs 
trait_fm = fread(paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",traitName,"_peaks_overlap.bed"),data.table = F,stringsAsFactors = F)
trait_fm.sub = trait_fm[,c("V1","V3","V4","V5","V9")]
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")
trait_fm.sub = subset(trait_fm,V5 > 0.2)[,c("V1","V3","V4","V5","V9")]
dim(trait_fm.sub)
nrow(subset(trait_fm,V5 > 0.9)[,c("V1","V3","V4","V5","V9")])
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")

# b: subset to those SNPs that lie in SCENT peaks
trait_fm = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/rbc_peaks_overlap.bed",data.table = F,stringsAsFactors = F)
trait_fm.sub = trait_fm[,c("V1","V3","V4","V5","V9")]
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")
trait_fm.sub = subset(trait_fm,V5 > 0.2)[,c("V1","V3","V4","V5","V9")]
dim(trait_fm.sub)
nrow(subset(trait_fm,V5 > 0.9)[,c("V1","V3","V4","V5","V9")])
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")

# 2. find the subset of peaks that are linked to genes
# df.sub = merge(subset(res.df.mg,!is.na(beta.y) & fdr.y < 0.1),subset(trait_fm.sub,pip>0.2),by="peak")
df.sub = merge(res.df.mg,trait_fm.sub,by="peak",all.x=TRUE)
df.sub = df.sub[,c("chr","snp_pos","rsid","pip","peak","gene","beta.x","pval.x","fdr.x","beta.y","pval.y","fdr.y")]
colnames(df.sub)[7:12] = c("beta_h","p_h","fdr_h","beta_t21","p_t21","fdr_t21")

# 3. 
small_ATAC_pb = fread("~/Downloads/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)
large_RNA_pb = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
small_RNA_pb = fread("~/Downloads/pb_de.hsc.txt",data.table = F,stringsAsFactors = F)
cyc_HSC_pb = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/HSC_v_CyclingHSC.txt",data.table = F,stringsAsFactors = F)
cyc_HSC_pb$P.Value.rank = rank(cyc_HSC_pb$P.Value)/nrow(cyc_HSC_pb)
colnames(small_ATAC_pb) = paste0(colnames(small_ATAC_pb),"_smATAC")
colnames(large_RNA_pb) = paste0(colnames(large_RNA_pb),"_lgRNA")
colnames(small_RNA_pb) = paste0(colnames(small_RNA_pb),"_smRNA")
colnames(cyc_HSC_pb) = paste0(colnames(cyc_HSC_pb),"_cycHSC")
df.sub2 = merge(df.sub,small_ATAC_pb,by.x='peak',by.y='names_smATAC',all.x=T)
df.sub2 = merge(df.sub2,small_RNA_pb,by.x='gene',by.y='names_smRNA',all.x=T)
df.sub2 = merge(df.sub2,large_RNA_pb,by.x='gene',by.y='names_lgRNA',all.x=T)
df.sub2 = merge(df.sub2,cyc_HSC_pb,by.x='gene',by.y='names_cycHSC',all.x=T)
df.sub2$fm = !(df.sub2$pip <= 0.2 | is.na(df.sub2$pip))

df.sub2[1,]
subset(df.sub2,fdr_t21 < 0.1 & adj.P.Val_smATAC < 0.1 & adj.P.Val_lgRNA < 0.1 & adj.P.Val_smRNA < 0.1 & fm)
subset(df.sub2,fdr_t21 < 0.1 & adj.P.Val_smATAC < 0.1 & adj.P.Val_lgRNA < 0.1 & P.Value_smRNA < 0.05 & fm)
subset(df.sub2,fdr_t21 < 0.2 & P.Value_smATAC < 0.1 & P.Value_lgRNA < 1 & P.Value_smRNA < 0.1 & fm=="TRUE")
subset(df.sub2,fdr_t21 < 0.1 & fm=="TRUE")
subset(df.sub2,fdr_t21 < 0.2 & fm=="TRUE" & P.Value_lgRNA < 0.05)
subset(df.sub2,fdr_t21 < 0.2 & fm=="TRUE" & P.Value_lgRNA < 0.2)

subset(df.sub2,fdr_t21 < 0.2 & fm=="TRUE" & adj.P.Val_cycHSC < 0.1)
subset(df.sub2,fdr_t21 < 0.2 & fm=="TRUE")

subset(df.sub2,gene=="ABCA7")

subset(df.sub2,gene=="TFR2")
subset(df.sub2,gene=="TFR2")[2,c("fdr_t21","adj.P.Val_smATAC","adj.P.Val_lgRNA","P.Value_smRNA","fm")]
subset(df.sub2,gene=="ABCA7")[,c("fdr_t21","adj.P.Val_smATAC","adj.P.Val_lgRNA","P.Value_smRNA","fm")]
subset(df.sub2,gene=="ABCA7")[,c("fdr_t21","P.Value_smATAC","P.Value_lgRNA","P.Value_smRNA","fm")]


# ggplot(df.sub2,aes(x=fm,y=-log10(p_t21))) + geom_boxplot()
# res = list()
# res[[1]] = fisher.test(df.sub2$p_t21 < 0.05,df.sub2$fm,alternative = 'greater')
# res[[2]] = fisher.test(df.sub2$fdr_t21 < 0.1,df.sub2$fm,alternative = 'greater')
# res[[3]] = fisher.test(df.sub2$fdr_t21 < 0.01,df.sub2$fm,alternative = 'greater')
# df.sub3 = subset(df.sub2,p_t21 < 0.05)
# res[[4]] = fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm,alternative = 'greater')
# table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
# df.sub3 = subset(df.sub2,fdr_t21 < 0.1)
# res[[5]] = fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm,alternative = 'greater')
# table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
# df.sub3 = subset(df.sub2,fdr_t21 < 0.01)
# res[[6]] = fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm,alternative = 'greater')

library("exact2x2")
res = list()
res[[1]] = exact2x2(df.sub2$p_t21 < 0.05,df.sub2$fm,alternative = 'two.sided')
res[[2]] = exact2x2(df.sub2$fdr_t21 < 0.2,df.sub2$fm,alternative = 'two.sided')
res[[3]] = exact2x2(df.sub2$fdr_t21 < 0.1,df.sub2$fm,alternative = 'two.sided')
# res[[3]] = exact2x2(df.sub2$fdr_t21 < 0.01,df.sub2$fm,alternative = 'two.sided')
df.sub3 = subset(df.sub2,p_t21 < 0.05)
res[[4]] = exact2x2(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm,alternative = 'two.sided')
table(df.sub3$P.Value_smRNA < 0.05,df.sub3$fm)
df.sub3 = subset(df.sub2,fdr_t21 < 0.2)
res[[5]] = exact2x2(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm,alternative = 'two.sided')
table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
subset(df.sub3,adj.P.Val_lgRNA < 0.1 & fm)
df.sub3 = subset(df.sub2,fdr_t21 < 0.1)
res[[6]] = exact2x2(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm,alternative = 'two.sided')
# table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
# df.sub3 = subset(df.sub2,fdr_t21 < 0.01)
# res[[6]] = exact2x2(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm,alternative = 'two.sided')

table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
# fisher.test(df.sub3$P.Value_lgRNA < 0.05 | df.sub3$P.Value_smRNA < 0.05,df.sub3$fm)
# fisher.test(df.sub3$P.Value_lgRNA < 0.05 | df.sub3$P.Value_smRNA < 0.05,df.sub3$fm)
df.sub3 = subset(df.sub2,fdr_t21 < 0.1)
table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
binom.test(2,5,p=66/(693+66))
binom.test(24,(24+67),p=6621/(6621+31320))
# density plots
# ggplot(df.sub3,aes(x=-log10(adj.P.Val_lgRNA),fill=df.sub3$fm)) + geom_density(adjust=1/2,alpha=0.4)
# ggplot(df.sub3,aes(x=(adj.P.Val_lgRNA),fill=df.sub3$fm)) + geom_density(adjust=0.5,alpha=0.4) + theme_bw() + theme(panel.grid = element_blank())
# ggplot(df.sub2,aes(x=(p_t21),fill=fm)) + geom_density(adjust=0.5,alpha=0.4) + theme_bw() + theme(panel.grid = element_blank())
# ggplot(df.sub2,aes(x=(fdr_t21),fill=fm)) + geom_density(adjust=0.5,alpha=0.4) + theme_bw() + theme(panel.grid = element_blank())
# ggplot(df.sub3,aes(x=(P.Value_lgRNA),fill=df.sub3$fm)) + geom_density(adjust=1/2,alpha=0.4)
# ggplot(df.sub3,aes(x=-log10(P.Value_lgRNA),fill=df.sub3$fm)) + geom_density(adjust=1/2,alpha=0.4)

# ggplot(df.sub3,aes(x=(P.Value_lgRNA),col=df.sub3$fm)) + geom_histogram(stat='frequency')
# ggplot(df.sub3,aes(x=(P.Value_lgRNA),col=df.sub3$fm)) + geom_histogram(aes(y = after_stat(count / sum(count))))

df.sub3 = subset(df.sub2,fdr_t21 < 0.1)
fisher.test(df.sub3$P.Value_lgRNA < 0.05 | df.sub3$P.Value_smATAC < 0.05,df.sub3$fm)
table(df.sub3$P.Value_lgRNA < 0.05 | df.sub3$P.Value_smATAC < 0.05,df.sub3$fm)

fisher.test(df.sub3$P.Value_lgRNA < 0.05 & df.sub3$P.Value_smATAC < 0.05,df.sub3$fm)
fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1 & df.sub3$adj.P.Val_smATAC < 0.1,df.sub3$fm)
fisher.test(df.sub3$adj.P.Val_lgRNA < 0.01 & df.sub3$adj.P.Val_smATAC < 0.01,df.sub3$fm)

estimate=unlist(lapply(res,function(x) x$estimate))
lower=unlist(lapply(res,function(x) x$conf.int[1]))
upper=unlist(lapply(res,function(x) x$conf.int[2]))
pvalue=unlist(lapply(res,function(x) x$p.value))
res <- data.frame(estimate, lower, upper,pvalue,row.names = NULL)
res$analy = NA
res$analy[1:3] = "Active cis-regulatory region\n(in T21 HSCs)"
res$analy[4:6] = "Differential expression of target GWAS genes\n(in large scRNA-seq HSCs)"
# res$analy[1:3] = "SNPs in active cis-regulatory region\n(in T21 HSCs)"
# res$analy[4:6] = "Targets genes are differentially expressed in T21\n(in large scRNA-seq HSCs)"
res$thres = NA
res$thres[c(1,4)] = "P < 0.05"
res$thres[c(2,5)] = "FDR < 0.2"
res$thres[c(3,6)] = "FDR < 0.1"
res$analy = factor(res$analy,levels=rev(unique(res$analy)))
# levels(res$analy) = rev(unique(res$analy))
library(ggplot2)

g=ggplot(subset(res,thres!="FDR < 0.01"),aes(x=analy,y=estimate,ymin=lower,ymax=upper,col=thres)) + 
  scale_color_manual(values=c("Blue","Purple","DarkRed")) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=90)) + 
  coord_flip() +
  geom_errorbar(width=0.3,position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5),size=rel(3)) +
  labs(x="",y="RBC GWAS Enrichment (OR)",col="Peak-gene threshold") +
  geom_hline(yintercept = 1,col='red',lty='dashed',lwd=1) +
  scale_y_continuous(trans="log10");g

f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.",traitName,".SCENT_enrich.boot.pdf")
pdf(f.plot,height=1.37*2,width=4.2*2)
print(g)
dev.off()

# binom.test(4,16,p=225/(3630+225))

