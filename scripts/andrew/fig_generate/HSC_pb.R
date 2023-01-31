library(ggplot2)
library(edgeR)
library(data.table)
geneName="TFR2"
# out = readRDS("~/Documents/Research/t21-proj/out/full/data_pb_leiden/Liver.pb.HSCs_MPPs.txt")
# df.aggre = out[[1]]
# metadata_to_use = out[[2]]
df.aggre = fread("~/Documents/Research/t21-proj/out/full/data_pb_leiden/Liver.pb.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
rownames(df.aggre) =df.aggre[,1]; df.aggre =  df.aggre[,-1]
df.aggre.cpm = cpm(df.aggre)

out = data.frame(disease=metadata_to_use$environment,expr=as.numeric(df.aggre.cpm[geneName,]))
out$disease = ifelse(out$disease=="DownSyndrome","T21","H")
library(ggplot2)
col_to_use = c("#FB8072" ,"#CBD5E8")
g=ggplot(out,aes(x=disease,y=expr,fill=disease)) +
  theme_bw() +
  # geom_violin() +
  geom_boxplot(outlier.shape=NA,width=0.2) +
  geom_point(alpha=0.6) +
  scale_fill_manual(values=col_to_use) +
  # scale_fill_brewer(palette = "Set1") +
  theme(panel.grid = element_blank()) + 
  labs(x="Disease",y=bquote(italic(.(geneName))~"expression in large scRNA-seq (cpm)")) + coord_flip() +
  guides(fill="none") +
  scale_y_continuous(trans = "log10") +
  ggpubr::theme_pubr(); g

f.plot = paste0("~/Documents/Research/t21_multiome/output/scent/plots/HSC_pb.",geneName,".pdf")
pdf(f.plot,height=2.5,width=5)
print(g)
dev.off()

g = g + geom_violin() + geom_boxplot(outlier.shape=NA,width=0.2) +   geom_point(alpha=0.6)
f.plot = paste0("~/Documents/Research/t21_multiome/output/scent/plots/HSC_pb.",geneName,".violin.pdf")
# pdf(f.plot,height=8,width=10)
pdf(f.plot,height=3,width=5)
print(g)
dev.off()

df.aggre2 = fread("~/Documents/Research/t21-proj/out/full/data_pb_leiden/Liver.pb.Cycling HSCs_MPPs.txt",data.table = F,stringsAsFactors = F)
rownames(df.aggre2) = df.aggre2[,1]
df.aggre2 = df.aggre2[,-1]
df.aggre2.cpm = cpm(df.aggre2)
out2 = rbind(out,
             data.frame(disease="T21 Cycling",expr=as.numeric(df.aggre.cpm[geneName,])))
g=ggplot(out2,aes(x=disease,y=expr,fill=disease)) +
  theme_bw() +
  # geom_violin() +
  geom_boxplot(outlier.shape=NA,width=0.2) +
  geom_point(alpha=0.6) +
  # scale_fill_manual(values=col_to_use) +
  # scale_fill_brewer(palette = "Set1") +
  theme(panel.grid = element_blank()) + 
  labs(x="Disease",y=bquote(italic(.(geneName))~"expression in large scRNA-seq (cpm)")) + coord_flip() +
  guides(fill="none") +
  scale_y_continuous(trans = "log10") +
  ggpubr::theme_pubr(); g
plot(1:5)



