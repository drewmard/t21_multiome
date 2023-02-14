library("exact2x2")
library(data.table)

df.sub2 = fread("~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs.DE.int.txt",data.table = F,stringsAsFactors = F)

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

traitName="rbc"
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.",traitName,".SCENT_enrich.boot.pdf")
pdf(f.plot,height=1.37*2,width=4.2*2)
print(g)
dev.off()
