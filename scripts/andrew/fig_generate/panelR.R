# 0. SCENT!
celltype_to_use="HSCs_H"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
res.df.h = fread(output_file,data.table = F,stringsAsFactors = F)

celltype_to_use="HSCs_T21"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
res.df.t21 = fread(output_file,data.table = F,stringsAsFactors = F)

res.df.mg = merge(res.df.h,res.df.t21,all=TRUE,by=c("gene","peak"))

# 1. find trait-associated SNPs 
trait_fm = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/rbc_peaks_overlap.bed",data.table = F,stringsAsFactors = F)
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
df.sub = df.sub[,c("chr","snp_pos","rsid","pip","peak","gene","beta.x","p.x","fdr.x","beta.y","p.y","fdr.y")]
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

df.sub2[order(df.sub2$pip,decreasing = T),]

ggplot(df.sub2,aes(x=fm,y=-log10(p_t21))) + geom_boxplot()
res = list()
res[[1]] = fisher.test(df.sub2$p_t21 < 0.05,df.sub2$fm)
res[[2]] = fisher.test(df.sub2$fdr_t21 < 0.1,df.sub2$fm)
res[[3]] = fisher.test(df.sub2$fdr_t21 < 0.01,df.sub2$fm)

df.sub3 = subset(df.sub2,p_t21 < 0.05)
res[[4]] = fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
df.sub3 = subset(df.sub2,fdr_t21 < 0.1)
res[[5]] = fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
df.sub3 = subset(df.sub2,fdr_t21 < 0.01)
res[[6]] = fisher.test(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
# fisher.test(df.sub3$P.Value_lgRNA < 0.05 | df.sub3$P.Value_smRNA < 0.05,df.sub3$fm)
# fisher.test(df.sub3$P.Value_lgRNA < 0.05 | df.sub3$P.Value_smRNA < 0.05,df.sub3$fm)
table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)
binom.test(4,16,p=225/(3630+225))

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
res$thres[c(2,5)] = "FDR < 0.1"
res$thres[c(3,6)] = "FDR < 0.01"
res$analy = factor(res$analy,levels=rev(unique(res$analy)))
# levels(res$analy) = rev(unique(res$analy))
library(ggplot2)
# ggplot(res,aes(x=analy,y=estimate,ymin=lower,ymax=upper,fill=thres)) + 
#   geom_bar(stat='identity',col='black',position=position_dodge()) + 
#   scale_fill_brewer(palette="Pastel1") +
#   ggpubr::theme_pubr() +
#   theme(axis.text.x = element_text(angle=90)) + 
#   coord_flip() +
#   geom_errorbar(width=0.2,position=position_dodge(width=0.9)) +
#   geom_point(position=position_dodge(width=0.9)) +
#   labs(x="",y="RBC GWAS Enrichment (OR)",fill="Peak-gene threshold") +
#   geom_hline(yintercept = 1,col='red',lty='dashed',lwd=1) +
#   scale_y_continuous(trans="log10")
  
g=ggplot(res,aes(x=analy,y=estimate,ymin=lower,ymax=upper,col=thres)) + 
  scale_color_manual(values=c("Blue","Purple","DarkRed")) +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(angle=90)) + 
  coord_flip() +
  geom_errorbar(width=0.2,position=position_dodge(width=0.9)) +
  geom_point(position=position_dodge(width=0.9)) +
  labs(x="",y="RBC GWAS Enrichment (OR)",col="Peak-gene threshold") +
  geom_hline(yintercept = 1,col='red',lty='dashed',lwd=1) +
  scale_y_continuous(trans="log10");g

traitName="rbc"
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.",traitName,".SCENT_enrich.pdf")
pdf(f.plot,height=1.37*2,width=4.2*2)
print(g)
dev.off()



subset(df.sub3,fm==1 & gene=="HBS1L")
subset(df.sub3,fm==1 & adj.P.Val_lgRNA < 0.1 & adj.P.Val_smATAC < 0.1)

table(df.sub3$adj.P.Val_lgRNA < 0.1,df.sub3$fm)

fisher.test(df.sub2$adj.P.Val_lgRNA < 0.1,df.sub2$fm)

fisher.test(df.sub2$P.Value_smATAC < 0.05,df.sub2$fm)
fisher.test(df.sub2$adj.P.Val_smATAC < 0.05,df.sub2$fm)
fisher.test(df.sub2$adj.P.Val_smRNA < 0.05,df.sub2$fm)
fisher.test(df.sub2$adj.P.Val_cycHSC < 0.1,df.sub2$fm)
table(df.sub2$P.Value.rank_cycHSC < 0.1,df.sub2$fm)

rank(df.sub2$P.Value_cycHSC)/nrow(df.sub2)
table(df.sub2$adj.P.Val_cycHSC < 0.0001)

table(df.sub2$adj.P.Val_lgRNA < 0.1,df.sub2$fm)
fisher.test(df.sub2$adj.P.Val_lgRNA < 0.1,df.sub2$fm)

table(df.sub2$fdr_t21 < 0.05,df.sub2$fm)

# fisher.test(df.sub2$p_h < 0.05,df.sub2$fm)

gene_of_interest = "TFR2"; peak_of_interest = "chr7-100639909-100642992"
gene_of_interest = "ABCA7"; peak_of_interest = "chr19-1063288-1063650"
Idents(dfseurat) <- "kmeans_RNA"
levels(dfseurat) = c("1","2","3")
DefaultAssay(dfseurat) = "ATAC"
col_to_use = c("#FB8072" ,"#CBD5E8","violetred3")
g = CoveragePlot(
  object = dfseurat,
  region = peak_of_interest,
  features = gene_of_interest,
  extend.upstream = 100,
  extend.downstream = 100
) & scale_fill_manual(values=col_to_use)
f.plot = paste0("/home/amarder/tmp/peak_coverage_plot.",gene_of_interest,".pdf")
pdf(f.plot,height=4,width=10)
print(g)
dev.off()

dfseurat[]

DefaultAssay(dfseurat) = "RNA"
FindMarkers(dfseurat,features=c("TFR2","ABCA7"),ident.1="2",ident.2="1")
FindMarkers(dfseurat,features=c("TFR2","ABCA7"),ident.1="2",ident.2=c("1","3"))

# link peaks to genes
y <- LinkPeaks(
  object = dfseurat[,dfseurat$disease=="T21"],
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = c("TFR2")
)


FindMarkers(dfseurat,features="TFR2",ident.1="2")



