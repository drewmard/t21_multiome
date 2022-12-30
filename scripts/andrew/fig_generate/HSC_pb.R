out = readRDS("~/Documents/Research/t21-proj/out/full/data_pb_leiden/Liver.pb.HSCs_MPPs.txt")
df.aggre = out[[1]]
metadata_to_use = out[[2]]
df.aggre.cpm = cpm(df.aggre)
out = data.frame(disease=metadata_to_use$environment,expr=as.numeric(df.aggre.cpm["TFR2",]))
out$disease = ifelse(out$disease=="DownSyndrome","T21","H")
library(ggplot2)
g=ggplot(out,aes(x=disease,y=expr,fill=disease)) +
  theme_bw() +
  geom_boxplot(outlier.shape=NA) +
  geom_point(alpha=0.6) +
  scale_fill_brewer(palette = "Set1") +
  theme(panel.grid = element_blank()) + 
  labs(x="Disease",y=bquote("HSC pseudobulk expression (log10 cpm) in large scRNA-seq")) + coord_flip() +
  guides(fill="none") +
  scale_y_continuous(trans = "log10") +
  ggpubr::theme_pubr(); g


# colnames(df.aggre.cpm)==rownames(metadata_to_use)
