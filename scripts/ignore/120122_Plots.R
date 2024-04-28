traitName="rbc"

library(RColorBrewer)
set.seed(031995)
n <- length(unique(trait_mat.save$subclust_v6))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_to_use = sample(col_vector, n)

viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF",
            "#5DC863FF", "#AADC32FF", "#FDE725FF")

res=list()
for (traitName in c("rbc","wbc","lymph")) {
  
  # traitName="rbc"
  trait_mat.save = fread(paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt"),data.table = F,stringsAsFactors = F)

  print(traitName)
  hsc.sub = subset(trait_mat.save,subclust_v6=="HSCs")
  # median(hsc.sub$TRS[hsc.sub$disease=="T21"])
  # median(hsc.sub$TRS[hsc.sub$disease=="H"])
  # aggregate(hsc.sub$TRS,by=list(hsc.sub$disease),median)
  # tmp=aggregate(hsc.sub$TRS,by=list(hsc.sub$dataset,hsc.sub$disease),median)
  # wilcox.test(subset(tmp,Group.2=="H")$x,subset(tmp,Group.2=="T21")$x)
  # tmp=aggregate(hsc.sub$TRS,by=list(hsc.sub$dataset,hsc.sub$disease),var)
  # tmp=aggregate(hsc.sub$TRS,by=list(hsc.sub$dataset,hsc.sub$disease),function(x) mean(x > quantile(trait_mat.save$TRS,probs=0.8)))
  # print(wilcox.test(subset(tmp,Group.2=="H")$x,subset(tmp,Group.2=="T21")$x)$p.value)

  
  # tmp=aggregate(hsc.sub$TRS,by=list(hsc.sub$dataset,hsc.sub$disease),var)
  # print(wilcox.test(subset(tmp,Group.2=="H")$x,subset(tmp,Group.2=="T21")$x)$p.value)
  # # ggplot(tmp,aes(x=Group.2,y=x)) + geom_boxplot(outlier.shape=NA) + geom_point() + theme_bw() + labs(x="Disease",y="Within-sample TRS variance")
  # g=ggplot(tmp,aes(x=Group.2,y=x,fill=Group.2)) +
  #   theme_bw() +
  #   geom_boxplot(outlier.shape=NA) +
  #   geom_point(alpha=0.6) +
  #   scale_fill_brewer(palette = "Set1") +
  #   theme(panel.grid = element_blank()) +
  #   labs(x="Disease",y=paste0("SCAVENGE TRS (",traitName,") variance within sample")) + coord_flip() +
  #   guides(fill="none")
  # print(g)
  # 
  # f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".HSCs.variance.pdf")
  # pdf(f.plot,width = 5,height=1.5)
  # print(g)
  # dev.off()
  
  # tmp=aggregate(hsc.sub$TRS,by=list(hsc.sub$dataset,hsc.sub$disease),median)
  # print(wilcox.test(subset(tmp,Group.2=="H")$x,subset(tmp,Group.2=="T21")$x)$p.value)
  # # ggplot(tmp,aes(x=Group.2,y=x)) + geom_boxplot(outlier.shape=NA) + geom_point() + theme_bw() + labs(x="Disease",y="Within-sample TRS variance")
  # g=ggplot(tmp,aes(x=Group.2,y=x,fill=Group.2)) +
  #   theme_bw() +
  #   geom_boxplot(outlier.shape=NA) +
  #   geom_point(alpha=0.6) +
  #   scale_fill_brewer(palette = "Set1") +
  #   theme(panel.grid = element_blank()) +
  #   labs(x="Disease",y=paste0("SCAVENGE TRS (",traitName,") median within sample")) + coord_flip() +
  #   guides(fill="none")
  # print(g)
  # 
  # f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".HSCs.med.pdf")
  # pdf(f.plot,width = 5,height=1.5)
  # print(g)
  # dev.off()
  
  
  
  
  # wilcox.test(tmp$x[1:6],tmp$x[7:12])
  # wilcox.test(hsc.sub$TRS[hsc.sub$disease=="H"],hsc.sub$TRS[hsc.sub$disease=="T21"])
  # 
  # wilcox.test(hsc.sub$TRS[hsc.sub$disease=="H"],hsc.sub$TRS[hsc.sub$disease=="T21"])
  print(res[[traitName]] <- fisher.test(hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9),
                    hsc.sub$disease))
  
  print(table(hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9),
                    hsc.sub$disease))
  
  
  
  # hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9)
  # print(aggregate(hsc.sub$TRS > quantile(trait_mat.save$TRS,probs=0.9),by=list(disease=hsc.sub$disease),mean))
  # 
  
  # # TRS boxplot:
  # g=ggplot(trait_mat.save,aes(x=subclust_v6,y=TRS,fill=subclust_v6)) +
  # theme_bw() +
  # geom_boxplot() +
  # scale_fill_brewer(palette = "Set1") +
  # theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  # labs(x="Cell type",y=paste0("TRS (",traitName,")")) +
  # scale_fill_manual(values=col_to_use) +
  # guides(fill=F)
  # print(g)
  # 
  # f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".no_disease.pdf")
  # pdf(f.plot,width = 10,height=4)
  # print(g)
  # dev.off()

  # g=ggplot(trait_mat.save,aes(x=subclust_v6,y=TRS,fill=subclust_v6)) +
  #   theme_bw() +
  #   geom_boxplot(outlier.shape = NA) +
  #   scale_fill_brewer(palette = "Set1") +
  #   theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
  #   labs(x="Cell type",y=paste0("TRS (",traitName,")")) +
  #   scale_fill_manual(values=col_to_use) +
  #   guides(fill=F)
  # # print(g)
  # 
  # f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".no_disease_no_pts.pdf")
  # pdf(f.plot,width = 10,height=4)
  # print(g)
  # dev.off()
  
  # # TRS U-MAP:
  # p2 <- ggplot(data=trait_mat.save, aes(UMAP_1, UMAP_2, color=TRS)) + 
  #   geom_point(size=1, na.rm = TRUE, alpha = 0.6) +
  #   scale_color_gradientn(colors = viridis) + 
  #   scale_alpha()+
  #   theme_bw() + 
  #   theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
  # # p2
  # f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/",traitName,".2.pdf")
  # pdf(f.plot,width=5,height=4)
  # print(p2)
  # dev.off()
  
  # # HSC boxplot:
  # g=ggplot(subset(trait_mat.save,subclust_v6=="HSCs"),aes(x=disease,y=TRS,fill=disease)) +
  #   theme_bw() +
  #   geom_boxplot() +
  #   scale_fill_brewer(palette = "Set1") +
  #   theme(panel.grid = element_blank()) +
  #   labs(x="Disease",y=paste0("SCAVENGE TRS (",traitName,")")) + coord_flip() +
  #   guides(fill="none")
  # print(g)
  # 
  # f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".HSCs.pdf")
  # pdf(f.plot,width = 5,height=1.5)
  # print(g)
  # dev.off()
  
}

res$lymph$conf.int[1]
estimate=unlist(lapply(res,function(x) x$estimate))
lower=unlist(lapply(res,function(x) x$conf.int[1]))
upper=unlist(lapply(res,function(x) x$conf.int[2]))
pvalue=unlist(lapply(res,function(x) x$p.value))
df <- data.frame(celltype=names(res), estimate, lower, upper,pvalue,row.names = NULL)

library(ggplot2)
fp <- ggplot(data=df, aes(x=celltype, y=estimate, ymin=lower, ymax=upper)) +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cell type") + ylab("Odds (95% CI)") +
  theme_bw()  # use a white background

library(ggplot2)    
fp <- ggplot(data=df, aes(x=celltype, y=estimate, ymin=lower, ymax=upper,fill=celltype)) +
  geom_bar(stat = 'identity',col='black') +
  geom_pointrange() + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Cell type") + ylab("Odds Ratio (95% CI)") + ggtitle("Association of extreme GWAS-enriched HSCs (>90% TRS) and T21 status") +
  guides(fill=F) +
  ggpubr::theme_pubr() +
  scale_fill_manual(values=col_to_use[1:3]) +
  theme(plot.title = element_text(hjust=0.5))# +
  # coord_flip()
f.plot = paste0("~/Documents/Research/t21_multiome/output/scavenge_plots/gwas_enriched_cells.pdf")
pdf(f.plot,width = 5.16*1.5,height=2.7*1.5)
print(fp)
dev.off()



