peak.df = fread("~/Downloads/rbc.fm_peak.pip_0.9.txt",data.table = F,stringsAsFactors = F)
peak.df$accMax = apply(peak.df[,-1],1,max)
peak.df.sub = subset(peak.df,accMax > 0.05)
dim(peak.df.sub)
res1=binom.test(sum(peak.df$acc1 > 0.05),nrow(peak.df),mean(peak.df$acc3>0.05))
res2=binom.test(sum(peak.df$acc2 > 0.05),nrow(peak.df),mean(peak.df$acc3>0.05))
tmp = rbind(data.frame(k=1,p=res1$estimate,l=res1$conf.int[1],h=res1$conf.int[2]),
            data.frame(k=2,p=res2$estimate,l=res2$conf.int[1],h=res2$conf.int[2]))

col_to_use = c("#FB8072" ,"#CBD5E8")
ggplot(tmp,aes(x=as.factor(k),y=p*100,fill=as.factor(k))) +
  geom_bar(stat='identity',col='black',alpha=0.6) +
  geom_point() +
  geom_errorbar(aes(ymin=l*100,ymax=h*100),width=0.2) +
  geom_hline(yintercept=100*mean(peak.df$acc3>0.05),col='red',lty='dashed') +
  ggpubr::theme_pubr() +
  labs(x="Cluster",y="% of Fine-mapped Peaks Accessible (>5%)",title = "PIP > 0.9") +
  theme(plot.title=element_text(hjust=0.5)) +
  scale_fill_manual(values=col_to_use) +
  guides(fill='none') +
  ylim(0,100)


peak.df = fread("~/Downloads/rbc.fm_peak.pip_0.2.txt",data.table = F,stringsAsFactors = F)
peak.df$accMax = apply(peak.df[,-1],1,max)
peak.df.sub = subset(peak.df,accMax > 0.05)
dim(peak.df.sub)
res1=binom.test(sum(peak.df$acc1 > 0.05),nrow(peak.df),mean(peak.df$acc3>0.05))
res2=binom.test(sum(peak.df$acc2 > 0.05),nrow(peak.df),mean(peak.df$acc3>0.05))
tmp = rbind(data.frame(k=1,p=res1$estimate,l=res1$conf.int[1],h=res1$conf.int[2]),
            data.frame(k=2,p=res2$estimate,l=res2$conf.int[1],h=res2$conf.int[2]))

col_to_use = c("#FB8072" ,"#CBD5E8")
ggplot(tmp,aes(x=as.factor(k),y=p*100,fill=as.factor(k))) +
  geom_bar(stat='identity',col='black',alpha=0.6) +
  geom_point() +
  geom_errorbar(aes(ymin=l*100,ymax=h*100),width=0.2) +
  geom_hline(yintercept=100*mean(peak.df$acc3>0.05),col='red',lty='dashed') +
  ggpubr::theme_pubr() +
  labs(x="Cluster",y="% of Fine-mapped Peaks Accessible (>5%)",title = "PIP > 0.2") +
  theme(plot.title=element_text(hjust=0.5)) +
  scale_fill_manual(values=col_to_use) +
  guides(fill='none') +
  ylim(0,100)


#

# apply(peak.df[,2:4],2,median)
# peak.df.melt = melt(peak.df,by="peakid")
# ggplot(peak.df.melt,aes(x=variable,y=value)) + geom_boxplot()
# 
# celltype_to_use="HSCs_H"
# output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
# res.df.h = fread(output_file,data.table = F,stringsAsFactors = F)
# 
# celltype_to_use="HSCs_T21"
# output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
# res.df.t21 = fread(output_file,data.table = F,stringsAsFactors = F)
# 
# res.df.mg = merge(res.df.h,res.df.t21,all=TRUE,by=c("gene","peak"))
# 
# length(unique(subset(res.df.mg,(p.x < 0.05 | p.y < 0.05) & peak %in% peak.df$peak)$peak))
# length(unique(subset(res.df.mg,(p.x < 0.05 | p.y < 0.05))$peak))
# length(unique(res.df.mg$peak))
# 
# length(unique(res.df.mg$peak))
# 
# res1=binom.test(sum(peak.df$acc1 > 0.05),nrow(peak.df),mean(peak.df$acc3>0.05))
# 
# 
# 
# apply(peak.df[,2:4] > 0.05,2,mean)
# cor.test(peak.df$acc1,peak.df$acc2)
# cor.test(peak.df$acc1,peak.df$acc3)
# 
# 
# apply(peak.data.df[,2:4],2,median)
# apply(peak.data.df[,2:4] > 0.05,2,mean)
# 
# 
# 
# 
