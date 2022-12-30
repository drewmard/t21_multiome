res = list()
traitName_list = c("")

fileDir="~/Documents/Research/t21_multiome/output/scavenge/"
flst = list.files(fileDir,pattern = "txt",include.dirs = FALSE)
traitName_list = gsub("\\..*","",flst)

# subset(kmean,kmeans_RNA==2)

table(trait_mat.save$subclust_v6,trait_mat.save$disease)
# for (celltype in c("NK cells","B cells","Megakaryocytes","Pro-inflammatory macrophages")) {
for (celltype in c("Late erythroid","Early erythroid")) {
res=list()
for (traitName in traitName_list) {
  print(traitName)
  # load TRS data:
  trait_mat.save = fread(paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt"),data.table = F,stringsAsFactors = F)

  # data subset for panels H, I, + J:
  hsc.sub = subset(trait_mat.save,subclust_v6==celltype)
  # hsc.sub = subset(trait_mat.save,subclust_v6=="NK cells")
  # hsc.sub = subset(trait_mat.save,subclust_v6=="Megakaryocytes")
  # hsc.sub = subset(trait_mat.save,subclust_v6=="Pro-inflammatory macrophages")
  # hsc.sub = subset(trait_mat.save,subclust_v6=="B cells")
  
  (res[[traitName]] <- fisher.test(hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9),
                                   hsc.sub$disease))

}

estimate=unlist(lapply(res,function(x) x$estimate))
lower=unlist(lapply(res,function(x) x$conf.int[1]))
upper=unlist(lapply(res,function(x) x$conf.int[2]))
pvalue=unlist(lapply(res,function(x) x$p.value))
df <- data.frame(celltype=names(res), estimate, lower, upper,pvalue,row.names = NULL)
df

set.seed(031995)
n <- length(unique(trait_mat.save$subclust_v6))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, nrow(df))

fp <- ggplot(data=df, aes(x=celltype, y=estimate, ymin=lower, ymax=upper,fill=celltype)) +
  geom_bar(stat = 'identity',col='black') +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Cell type") + ylab("Odds Ratio (95% CI)") + ggtitle("Association of GWAS-enriched HSCs (>90% TRS) and T21 status") +
  guides(fill="none") +
  ggpubr::theme_pubr() +
  scale_fill_manual(values=col_to_use) +
  scale_y_continuous(trans="log10") +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle=30,hjust=1)); # print(fp)
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/gwas_enriched_cells.all.",celltype,".pdf")
pdf(f.plot,width = 9.1*1.3,height=3.2*1.3)
print(fp)
dev.off()
}

#########

kmean = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/data/hsc.kmeans.UMAP.txt",data.table = F,stringsAsFactors = F)
# subset(kmean,kmeans_RNA==2)

res=list()
for (traitName in traitName_list) {
  print(traitName)
  # load TRS data:
  trait_mat.save = fread(paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt"),data.table = F,stringsAsFactors = F)
  trait_mat.save = merge(trait_mat.save,kmean[,c('cell','dataset','kmeans_RNA')],by=c("cell","dataset"),all.x=TRUE)
  
  # data subset for panels H, I, + J:
  hsc.sub = subset(trait_mat.save,subclust_v6=="HSCs")
  
  # (res[[traitName]] <- fisher.test(hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9),
  #                                  hsc.sub$disease))
  (res[[traitName]] <- fisher.test(hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9),
                                   hsc.sub$kmeans_RNA==2))
  
}

estimate=unlist(lapply(res,function(x) x$estimate))
lower=unlist(lapply(res,function(x) x$conf.int[1]))
upper=unlist(lapply(res,function(x) x$conf.int[2]))
pvalue=unlist(lapply(res,function(x) x$p.value))
df <- data.frame(celltype=names(res), estimate, lower, upper,pvalue,row.names = NULL)
df

set.seed(031995)
n <- length(unique(trait_mat.save$subclust_v6))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, nrow(df))

fp <- ggplot(data=df, aes(x=celltype, y=estimate, ymin=lower, ymax=upper,fill=celltype)) +
  geom_bar(stat = 'identity',col='black') +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Cell type") + ylab("Odds Ratio (95% CI)") + ggtitle("Association of GWAS-enriched HSCs (>90% TRS) and cluster 2 HSCs") +
  guides(fill="none") +
  ggpubr::theme_pubr() +
  scale_fill_manual(values=col_to_use) +
  scale_y_continuous(trans="log10") +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x = element_text(angle=30,hjust=1)); fp
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/gwas_enriched_cells.c2.pdf")
pdf(f.plot,width = 9.1*1.3,height=3.2*1.3)
print(fp)
dev.off()

