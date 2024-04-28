library(data.table)
library(dplyr)
multiome_samples = c("15582 nuclei A","15582 nuclei B")
scrna_samples = c("T21 15582 B")

# multiome_samples = c("15669_A","15669_B")
# scrna_samples = c("15669H")

f="/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/Multiome.RNA_ATAC.meta.txt"
multiome = fread(f,data.table = F,stringsAsFactors = F)
annot = fread("~/Downloads/multiome_annot.csv",data.table = F,stringsAsFactors = F)
multiome = merge(multiome,annot,by="subclust_v6")
multiome.sub = subset(multiome,dataset %in% multiome_samples)
tab = aggregate(multiome.sub$subclust_v6,by=list(annot=multiome.sub$broad,disease=multiome.sub$disease,patient_sample=multiome.sub$dataset),length)
tab <- tab %>%
  group_by(disease,patient_sample) %>%
  mutate(prop = x / sum(x)) %>%
  as.data.frame()

f = "~/Documents/Research/t21-proj/out/full/data/meta.10X_DownSyndrome_Liver.umap2d.cells_removed.txt"
scrna = fread(f,data.table = F,stringsAsFactors = F)
scrna.sub = subset(scrna,sample==scrna_samples)

tab2 = aggregate(scrna.sub$leiden_names,by=list(annot=scrna.sub$cell_type_groups,patient_sample=scrna.sub$sample),length)
tab2 <- tab2 %>%
  group_by(patient_sample) %>%
  mutate(prop = x / sum(x)) %>%
  as.data.frame()

tab.mg = merge(tab2,tab,by="annot")

library(ggplot2)
g1 = ggplot(subset(tab.mg,patient_sample.y==multiome_samples[1]),aes(x=prop.x,y=prop.y,col=annot)) +
  geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') +
  scale_color_brewer(palette = 'Set2') +
  ggpubr::theme_pubr() +
  labs(x="scRNA proportion",y="Multiome proportion",title=paste0("Multiome: ",multiome_samples[1])) +
  labs(col="") +
  # guides(col="none") +
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)); g1

g2 = ggplot(subset(tab.mg,patient_sample.y==multiome_samples[2]),aes(x=prop.x,y=prop.y,col=annot)) +
  geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') +
  scale_color_brewer(palette = 'Set2') +
  ggpubr::theme_pubr() +
  labs(x="scRNA proportion",y="Multiome proportion",title=paste0("Multiome: ",multiome_samples[2])) +
  # guides(col="none") +
  labs(col="") +
  theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5)); g2

library(cowplot)
pdf("~/Downloads/15582.legend.pdf",width = 10*1.15,height=5*1.15)
# pdf("~/Downloads/15669.legend.pdf",width = 10*1.15,height=5*1.15)
print(plot_grid(g1,g2,ncol=2))
dev.off()
pdf("~/Downloads/15582.no_legend.pdf",width = 7,height=3)
# pdf("~/Downloads/15669.no_legend.pdf",width = 7,height=3)
print(plot_grid(g1+guides(col="none"),g2+guides(col="none"),ncol=2))
dev.off()

plot_grid(g1+guides(col="none"),g2+guides(col="none"),ncol=2)



cor(subset(tab.mg,patient_sample.y==multiome_samples[1])[,c("prop.x","prop.y")])
cor(subset(tab.mg,patient_sample.y==multiome_samples[2])[,c("prop.x","prop.y")])

cor.test(subset(tab.mg,patient_sample.y=="15582 nuclei A")$prop.x,
         subset(tab.mg,patient_sample.y=="15582 nuclei A")$prop.y)



