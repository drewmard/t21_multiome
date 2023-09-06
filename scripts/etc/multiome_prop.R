library(data.table)
f="/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/Multiome.RNA_ATAC.meta.txt"
df = fread(f,data.table = F,stringsAsFactors = F)
# fwrite(data.frame(unique(df$subclust_v6)),"~/Downloads/multiome_annot.csv",quote = F)
annot = fread("~/Downloads/multiome_annot.csv",data.table = F,stringsAsFactors = F)
df = merge(df,annot,by="subclust_v6")
tab2 = aggregate(df$subclust_v6,by=list(annot=df$broad,disease=df$disease,patient_sample=df$dataset),length)
library(dplyr)
tab2 <- tab2 %>%
  group_by(disease,patient_sample) %>%
  mutate(prop = x / sum(x))
tab2 = as.data.frame(tab2)
# tab2

library(ggplot2)
tab2$annot = factor(tab2$annot,levels=c("HSC/Progenitors", "Erythroid", "Mast cells", "Megakaryocytes", "B cells", "NK/T cells", "Myeloid", "Stroma"))
library(dplyr)
tab2$disease = recode(tab2$disease,"H" = "Disomic", "T21" = "Ts21")
g=ggplot(tab2,aes(x=disease,y=prop,fill=disease)) + 
  geom_violin() + geom_point() + 
  facet_grid(cols=vars(annot)) + ggpubr::theme_pubr()  + 
  scale_fill_brewer(palette="Set2") + guides(fill="none") +
  labs(y="Cell proportions",x="Trisomy status")
pdf("~/Downloads/multiome_prop.pdf",width=13,height=6)
print(g)
dev.off()



uniq = unique(annot$broad)
pval = rep(NA,length(uniq))
for (i in 1:length(uniq)) {
  cell=uniq[i]
  tab2.sub = subset(tab2,annot==cell)
  pval[i] = wilcox.test(tab2.sub$prop[tab2.sub$disease=="H"],tab2.sub$prop[tab2.sub$disease=="T21"])$p.value
}

res = data.frame(uniq,pval,fdr=p.adjust(pval,method='fdr'))
fwrite(res,"~/Downloads/multiome_celltype_compare.txt",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

library(ggplot2)
g2<-ggplot(tab2, aes(fill=annot, y=prop*100, x=reorder(patient_sample,disease=="T21"),col=disease)) + 
  geom_bar(position="stack", stat="identity") +
  scale_color_manual(values=c("orange","black"),labels=c("Disomic","Ts21"),name="Disease Status") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(hjust=1,angle=60),
        plot.title = element_text(hjust=0.5)) +
  scale_fill_brewer(palette="Set3") +
  labs(x="Patient ID",y="Cell composition (%)",title="CD45+",fill="Cell type");g2
