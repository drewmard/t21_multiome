library(data.table)
library(ggplot2)
library(RColorBrewer)

# load TRS data:
traitName="rbc"
trait_mat.save = fread(paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt"),data.table = F,stringsAsFactors = F)

# A
library(RColorBrewer)
library(ggtext)
tmp = trait_mat.save
set.seed(03191995)
n <- length(unique(tmp$subclust_v6))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, n)
col_to_use[9] = "#E5C494"
tmpout = aggregate(tmp[,c("UMAP_1","UMAP_2")],by=list(tmp$subclust_v6),mean)
p2 <- ggplot(data=tmp[sample(1:nrow(tmp),nrow(tmp),replace = F),], aes(UMAP_1, UMAP_2)) + 
  geom_point(aes(color=subclust_v6),size=1, na.rm = TRUE, alpha = 0.6) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank(),
        panel.grid = element_blank()) + 
  xlab("UMAP 1") + ylab("UMAP 2") +
  scale_color_manual(values=col_to_use)
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/umap.clusters.labels.pdf")
pdf(f.plot,width=10.2,height=7)
print(p2 +   geom_richtext(data = tmpout,aes(label=Group.1,x=UMAP_1,y=UMAP_2)))
dev.off()
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/umap.clusters.text.pdf")
pdf(f.plot,width=10.2,height=7)
print(p2 +   geom_text(data = tmpout,aes(label=Group.1,x=UMAP_1,y=UMAP_2)))
dev.off()
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/umap.clusters.none.pdf")
pdf(f.plot,width=10.2,height=7)
print(p2)
dev.off()



# geom_richtext(label=tmp$Group.1,x=tmp$UMAP_1,y=tmp$UMAP_2)
# scale_color_gradientn(colors = viridis) + 

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/umap.clusters.pdf")
pdf(f.plot,width=10.2,height=7)
print(p2)
dev.off()

library(RColorBrewer)
set.seed(031995)
n <- length(unique(trait_mat.save$subclust_v6))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, n)

viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF",
            "#5DC863FF", "#AADC32FF", "#FDE725FF")


res=list()
res.hsc_med = list()
res.hsc_var = list()
for (traitName in c("rbc","wbc","lymph")) {
  
  print(traitName)
  
  # load TRS data:
  trait_mat.save = fread(paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt"),data.table = F,stringsAsFactors = F)
  
  # data subset for panels H, I, + J:
  hsc.sub = subset(trait_mat.save,subclust_v6=="HSCs")
  
  # B, C, + D:
  p2 <- ggplot(data=trait_mat.save, aes(UMAP_1, UMAP_2, color=TRS)) +
    geom_point(size=1, na.rm = TRUE, alpha = 0.6) +
    scale_color_gradientn(colors = viridis) +
    scale_alpha()+
    theme_bw() +
    theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(panel.grid = element_blank())
  f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/",traitName,".2.pdf")
  pdf(f.plot,width=5,height=4)
  print(p2)
  dev.off()
  
  
  # E, F, + G:
  y=aggregate(trait_mat.save$TRS,by=list(trait_mat.save$subclust_v6),median)
  trait_mat.save.mg = merge(trait_mat.save,y,by.x="subclust_v6",by.y="Group.1")
  # g=ggplot(trait_mat.save.mg,aes(x=subclust_v6,y=scale(TRS),fill=subclust_v6,alpha=x)) +
  g=ggplot(trait_mat.save.mg,aes(x=subclust_v6,y=scale(TRS),fill=subclust_v6,alpha=ifelse(x > 1,1,0.2))) +
    theme_bw() +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1") +
    theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
    labs(x="Cell type",y=paste0("SCAVENGE TRS (",traitName,")")) +
    scale_fill_manual(values=col_to_use) +
    guides(fill="none",alpha="none"); g
  # g=ggplot(trait_mat.save,aes(x=subclust_v6,y=scale(TRS),fill=subclust_v6)) +
  #   theme_bw() +
  #   geom_boxplot() +
  #   scale_fill_brewer(palette = "Set1") +
  #   theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  #   labs(x="Cell type",y=paste0("SCAVENGE TRS (",traitName,")")) +
  #   scale_fill_manual(values=col_to_use) +
  #   guides(fill="none"); g
  f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".no_disease.pdf")
  pdf(f.plot,width = 10,height=4)
  print(g)
  dev.off()

  # H
  tmp=aggregate(hsc.sub$TRS,by=list(hsc.sub$dataset,hsc.sub$disease),median)
  res.hsc_med[[traitName]] = wilcox.test(subset(tmp,Group.2=="H")$x,subset(tmp,Group.2=="T21")$x)$p.value
  g=ggplot(tmp,aes(x=Group.2,y=x,fill=Group.2)) +
    theme_bw() +
    geom_boxplot(outlier.shape=NA) +
    geom_point(alpha=0.6) +
    scale_fill_brewer(palette = "Set1") +
    theme(panel.grid = element_blank()) +
    labs(x="Disease",y=paste0("SCAVENGE TRS (",traitName,") median within sample")) + coord_flip() +
    guides(fill="none") +
    ggpubr::theme_pubr()

  f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".HSCs.med.pdf")
  pdf(f.plot,width = 5,height=1.5)
  print(g)
  dev.off()
  
  # I
  tmp=aggregate(hsc.sub$TRS,by=list(hsc.sub$dataset,hsc.sub$disease),var)
  res.hsc_var[[traitName]] = (wilcox.test(subset(tmp,Group.2=="H")$x,subset(tmp,Group.2=="T21")$x)$p.value)
  g=ggplot(tmp,aes(x=Group.2,y=x,fill=Group.2)) +
    theme_bw() +
    geom_boxplot(outlier.shape=NA) +
    geom_point(alpha=0.6) +
    scale_fill_brewer(palette = "Set1") +
    theme(panel.grid = element_blank()) +
    labs(x="Disease",y=paste0("SCAVENGE TRS (",traitName,") variance within sample")) + coord_flip() +
    guides(fill="none") +
    ggpubr::theme_pubr()

  f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".HSCs.variance.pdf")
  pdf(f.plot,width = 5,height=1.5)
  print(g)
  dev.off()
  
  # J (1/2)
  (res[[traitName]] <- fisher.test(hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9),
                                        hsc.sub$disease))
  
}

# J (2/2)
estimate=unlist(lapply(res,function(x) x$estimate))
lower=unlist(lapply(res,function(x) x$conf.int[1]))
upper=unlist(lapply(res,function(x) x$conf.int[2]))
pvalue=unlist(lapply(res,function(x) x$p.value))
df <- data.frame(celltype=names(res), estimate, lower, upper,pvalue,row.names = NULL)

set.seed(031995)
n <- length(unique(trait_mat.save$subclust_v6))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, n)

fp <- ggplot(data=df, aes(x=celltype, y=estimate, ymin=lower, ymax=upper,fill=celltype)) +
  geom_bar(stat = 'identity',col='black') +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Cell type") + ylab("Odds Ratio (95% CI)") + ggtitle("Association of GWAS-enriched HSCs (>90% TRS) and T21 status") +
  guides(fill="none") +
  ggpubr::theme_pubr() +
  scale_fill_manual(values=col_to_use[1:3]) +
  theme(plot.title = element_text(hjust=0.5))
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/gwas_enriched_cells.pdf")
pdf(f.plot,width = 5.16*1.3,height=2.7*1.3)
print(fp)
dev.off()


# K
expr.df = fread("~/Downloads/expr.df.txt",data.table = F,stringsAsFactors = F)
gene_of_interest_lst = c("GATA1","ITGA2B","MKI67","PROM1","SPINK2")
expr.df.melt = reshape2::melt(expr.df,id.vars=c("cell","disease","trs_ex"))
expr.df.melt = subset(expr.df.melt,variable!="PROCR")
expr.df.melt$trs_ex = ifelse(expr.df.melt$trs_ex,"GWAS-enriched HSCs (rbc)","Other HSCs")
g = ggplot(expr.df.melt,aes(x=variable,y=value,fill=trs_ex)) + 
  geom_boxplot() +
  ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(face = "italic")) +
  labs(x="Gene",y="Log-normalized expression",fill="") +
  scale_fill_manual(values=c(viridis[length(viridis)],viridis[3]))
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/gwas_enriched_cells.rna_markers.pdf")
pdf(f.plot,height = 1.6*2,width=6.3*2)
print(g)
dev.off()

for (gene_of_interest in gene_of_interest_lst) {
  print(gene_of_interest)
  print(wilcox.test(expr.df[expr.df$trs_ex==0,gene_of_interest],expr.df[expr.df$trs_ex==1,gene_of_interest])$p.value)
}

# L-P
library(data.table)
df = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hsc.kmeans.UMAP.txt",data.table = F,stringsAsFactors = F)

# L
p2 <- ggplot(data=df[sample(1:nrow(df),nrow(df),replace = F),], aes(UMAP_1_RNA, UMAP_2_RNA)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6,aes(col=as.factor(disease)) )+
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  xlab("UMAP1") + ylab("UMAP2") + labs(col="disease status") +
  # scale_color_manual(values=col_to_use) +
  scale_color_brewer(palette = "Set1") +
  ggpubr::theme_pubr()
p2
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.rna_umap.disease.pdf")
pdf(f.plot,width=4.3,height=4.7)
print(p2)
dev.off()

library(RColorBrewer)
brewer.pal(n=4,"Set1")

ggplot(data=data.frame(x=rnorm(100),y=rnorm(100)),aes(x,y)) + geom_point(col="#E41A1C")
ggplot(data=data.frame(x=rnorm(100),y=rnorm(100)),aes(x,y)) + geom_point(col="#377EB8")

ggplot(data=data.frame(x=rnorm(100),y=rnorm(100))) + geom_point(col="#377EB8")

# M
# col_to_use = c("#FB8072" ,"#CBD5E8","violetred3")
col_to_use = c("#FB8072" ,"#CBD5E8","#CD3278")
# col_to_use = c("#FB8072" ,"violetred3","#CBD5E8")
p2 <- ggplot(data=df[sample(1:nrow(df),nrow(df),replace = F),], aes(UMAP_1_RNA, UMAP_2_RNA)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6,aes(col=as.factor(kmeans_RNA)) )+
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  xlab("UMAP1") + ylab("UMAP2") + labs(col="cluster") +
  scale_color_manual(values=col_to_use) +
  ggpubr::theme_pubr()
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.rna_umap.rna_kmeans.pdf")
pdf(f.plot,width=4.3,height=4.7)
print(p2)
dev.off()

# N
p2 <- ggplot(data=df[sample(1:nrow(df),nrow(df),replace = F),], aes(UMAP_1_ATAC, UMAP_2_ATAC)) + 
  geom_point(size=1, na.rm = TRUE, alpha = 0.6,aes(col=as.factor(disease)) )+
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  xlab("UMAP1") + ylab("UMAP2") + labs(col="disease status") +
  # scale_color_manual(values=col_to_use) +
  scale_color_brewer(palette = "Set1") +
  ggpubr::theme_pubr()
p2
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.atac_umap.disease.pdf")
pdf(f.plot,width=4.3,height=4.7)
print(p2)
dev.off()

# O
col_to_use = c("#FB8072" ,"#CBD5E8","violetred3")
# col_to_use = c("#FB8072" ,"violetred3","#CBD5E8")
p2 <- ggplot(data=df[sample(1:nrow(df),nrow(df),replace = F),], aes(UMAP_1_ATAC, UMAP_2_ATAC)) + 
  # geom_point(aes(color=as.factor(kmeans_RNA)),size=rel(1.5), na.rm = TRUE, alpha = 0.6) +
  geom_point(size=1, na.rm = TRUE, alpha = 0.6,aes(col=as.factor(kmeans_RNA)) )+
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  xlab("UMAP1") + ylab("UMAP2") + labs(col="cluster") +
  scale_color_manual(values=col_to_use) +
  ggpubr::theme_pubr()
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.atac_umap.rna_kmeans.pdf")
pdf(f.plot,width=4.3,height=4.7)
print(p2)
dev.off()

# P
traitName='rbc'
trait_mat.save = fread(paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt"),data.table = F,stringsAsFactors = F)
trait_mat.save$trs_ex = trait_mat.save$TRS > quantile(trait_mat.save$TRS,probs=0.9)
trait_mat.save$trs_ex = ifelse(trait_mat.save$trs_ex,"GWAS-enriched HSCs","Other HSCs")
# fwrite(trait_mat.save[,c("cell","dataset","TRS","trs_ex")],"~/Downloads/rbc_trs.txt",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
df.mg = merge(df,trait_mat.save[,c("cell","dataset","TRS","trs_ex")],by=c("cell","dataset"))
# trait_mat.save = subset(trait_mat.save,subclust_v6=="HSCs")
# df.mg = subset(df,subclust_v6=="HSCs")
# nrow(trait_mat.save); nrow(df.mg)
# # discrepancies because singleton HSCs are removed in SCAVENGE!

col_to_use = c("#FDE725FF","#440154FF")
# p2 <- ggplot(data=df.mg[sample(1:nrow(df.mg),nrow(df.mg),replace = F),], aes(UMAP_1_RNA, UMAP_2_RNA, color=trs_ex)) + 
p2 <- ggplot(data=df.mg[order(as.numeric(as.factor(df.mg$trs_ex)),decreasing = T),], aes(UMAP_1_RNA, UMAP_2_RNA, color=trs_ex)) +
# p2 <- ggplot(data=df.mg[order(as.numeric(as.factor(df.mg$trs_ex)),decreasing = T),], aes(UMAP_1_ATAC, UMAP_2_ATAC, color=trs_ex)) +
  geom_point(size=1, na.rm = TRUE, alpha = 1) +
  scale_alpha()+
  theme_bw() + 
  theme(legend.title = element_blank()) + 
  xlab("UMAP1") + ylab("UMAP2") + labs(col="") +
  scale_color_manual(values=col_to_use) +
  ggpubr::theme_pubr()
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.rna_umap.gwas_enriched.pdf")
pdf(f.plot,width=4.3,height=4.7)
print(p2)
dev.off()

tmp = aggregate(as.numeric(trs_ex=="GWAS-enriched HSCs")~as.factor(kmeans_RNA),df.mg,mean)
p2 = ggplot(df.mg,aes(x=as.factor(kmeans_RNA),y=scale(TRS))) + 
  geom_jitter(width=0.1,aes(col=trs_ex),alpha=0.6) + 
  geom_boxplot(outlier.shape=NA,alpha=0.8,fill='gray94',lwd=1,width=0.2,col="#FB8072") + 
  ggpubr::theme_pubr() +
  labs(x="Cluster",y="SCAVENGE TRS (rbc)") +
  scale_color_manual(values=c(viridis[length(viridis)],viridis[3])) + labs(col="") +
  scale_x_discrete(labels=paste0(tmp[,1],"\n (",round(tmp[,2]*100,1),"%)")); p2
f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.",traitName,".gwas_enriched.boxplot.pdf")
pdf(f.plot,width=3.8,height=5)
print(p2)
dev.off()


# load TRS data:


