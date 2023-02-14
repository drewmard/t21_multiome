library(data.table)
df = fread("~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs.DE.txt",data.table = F,stringsAsFactors = F)
subset(df,gene=="RAPGEF2")[3,]


df$chr = as.numeric(substring(gsub("-.*","",df$peak),4))
tmp = subset(df,fdr_t21 < 0.2)
tab = table(as.factor(tmp$chr))
tab = tab/sum(tab)
df1 = as.data.frame(tab)
df1$type = "Ts21"
tmp = subset(df,fdr_H < 0.2)
tab = table(as.factor(tmp$chr))
tab = tab/sum(tab)
df2 = as.data.frame(tab)
df2[,2] = -1*df2[,2]
df2$type = "H"
df3 = as.data.frame(rbind(df1,df2))
df3$type = factor(df3$type,levels = c("Ts21","H"))
library(ggplot2)
col_to_use = rev(c("#FB8072" ,"#CBD5E8"))
g <- ggplot(df3,
            aes(x=Var1,
                y=Freq,
                fill=type)) +
  geom_bar(stat='identity',position='identity',col='black') +
  theme_bw() + 
  theme(
    panel.grid=element_blank(),
    legend.title = element_blank()
  ) + 
  geom_abline(slope=0,intercept=0,col='black') +
  labs(y='Proportion',x='Chromosome',title = "Chromosomal distribution of significant peak-gene links") +
  scale_fill_manual(values=col_to_use) +
  scale_y_continuous(breaks=seq(-0.15,0.14,by=0.03),labels = abs(seq(-0.15,0.14,by=0.03))) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        plot.title = element_text(hjust=0.5));g

pdf("~/Documents/Research/t21_multiome/output/scent/plots/peak_gene_dist.pdf",width = 10,height=4)
print(g)
dev.off()


