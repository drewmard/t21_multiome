library(data.table)
traitName="rbc"
traitName="wbc"
traitName="lymph"

for (traitName in c("rbc","wbc","lymph")) { 
  f = paste0("~/Documents/Research/t21_multiome/output/scavenge/",traitName,".txt")
  df = fread(f,data.table = F,stringsAsFactors = F)
  df$trs_ex = df$TRS > quantile(df$TRS,probs=0.9)
  df$trs_ex = ifelse(df$trs_ex,"GWAS-enriched HSCs","Other HSCs")
  
  lineage = fread("~/Downloads/Lineage_Probabilities_per_Cell.csv",data.table = F,stringsAsFactors = F)
  mapping = fread("~/Downloads/dataset_metadata_to_cell.csv",data.table = F,stringsAsFactors = F,header=T)
  which(lineage$V1!=mapping$V1)
  lineage = cbind(mapping,lineage[,-1])
  lineage$cell = unlist(lapply(strsplit(lineage$V1,"_"),function(x) x[1]))
  df2 = subset(df,subclust_v6=="HSCs")
  
  df.mg = merge(df2,lineage,by=c("cell","dataset"),all=T)
  
  linNum=1
  df.mg$lin = apply(df.mg[,paste0("Lineage_",c(1:3),"_prob")],1,which.max)
  aggregate(TRS~lin,df.mg,median)
  aggregate(trs_ex=="GWAS-enriched HSCs"~lin,df.mg,mean)
  
  wilcox.test(df.mg$TRS[df.mg$lin==1],df.mg$TRS[df.mg$lin==2])
  wilcox.test(df.mg$TRS[df.mg$lin==3],df.mg$TRS[df.mg$lin==2])
  wilcox.test(df.mg$TRS[df.mg$lin==1],df.mg$TRS[df.mg$lin==3])
  
  viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF",
              "#5DC863FF", "#AADC32FF", "#FDE725FF")
  library(ggplot2)
  # tmp = aggregate(as.numeric(trs_ex=="GWAS-enriched HSCs")~as.factor(lin),df.mg,mean)
  tmp = aggregate((TRS)~as.factor(lin),df.mg,median)
  
  # p2 = ggplot(subset(df.mg,!is.na(TRS)),aes(x=as.factor(lin),y=scale(TRS))) + 
  #   geom_jitter(width=0.1,aes(col=trs_ex),alpha=0.6) + 
  #   geom_boxplot(outlier.shape=NA,alpha=0.8,fill='gray94',lwd=1,width=0.2,col="#FB8072") + 
  #   ggpubr::theme_pubr() +
  #   labs(x="HSC Branch",y=paste0("SCAVENGE TRS (",traitName,")")) +
  #   scale_color_manual(values=c(viridis[length(viridis)],viridis[3])) + labs(col="") +
  #   scale_x_discrete(labels=paste0(tmp[,1],"\n (",round(tmp[,2]*100,1),"%)")); p2
  p2 = ggplot(subset(df.mg,!is.na(TRS)),aes(x=as.factor(lin),y=(TRS))) + 
    # geom_jitter(width=0.1,aes(col=scale(TRS)),alpha=0.6) + 
    # geom_boxplot(outlier.shape=NA,alpha=0.8,fill='gray94',lwd=1,width=0.2,col="#FB8072") + 
    geom_jitter(width=0.15,aes(col=(TRS)),alpha=0.6) + 
    geom_boxplot(outlier.shape=NA,alpha=0.8,fill='gray94',lwd=1,width=0.3,col="#FB8072") + 
    ggpubr::theme_pubr() +
    labs(x="HSC Branch",y=paste0("SCAVENGE TRS (",traitName,")")) +
    guides(col="none") +
    scale_color_gradient(high=viridis[length(viridis)],low=viridis[3]) + labs(col="") +
    scale_x_discrete(labels=paste0(tmp[,1])); p2
  f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.",traitName,".gwas_enriched.lin.boxplot.pdf")
  pdf(f.plot,width=3.8,height=5)
  print(p2)
  dev.off()
}


diseaseType="H"
for (diseaseType in c("H","T21")) { 
  p2 = ggplot(subset(df.mg,!is.na(TRS) & disease==diseaseType),aes(x=as.factor(lin),y=(TRS))) + 
    geom_jitter(width=0.15,aes(col=(TRS)),alpha=0.6) + 
    geom_boxplot(outlier.shape=NA,alpha=0.8,fill='gray94',lwd=1,width=0.3,col="#FB8072") + 
    ggpubr::theme_pubr() +
    labs(x="HSC Branch",y=paste0("SCAVENGE TRS (",traitName,")")) +
    guides(col="none") +
    scale_color_gradient(high=viridis[length(viridis)],low=viridis[3]) + labs(col="") +
    scale_x_discrete(labels=paste0(tmp[,1])); p2
  f.plot = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scavenge_plots/HSCs.",traitName,".gwas_enriched.lin.disease_",diseaseType,".boxplot.pdf")
  pdf(f.plot,width=3.8,height=5)
  print(p2)
  dev.off()
}


cor(df.mg$TRS,df.mg$Lineage_3_prob,use='na.or.complete')
cor(df.mg$TRS,df.mg$Lineage_2_prob,use='na.or.complete')
cor(df.mg$TRS,df.mg$Lineage_1_prob,use='na.or.complete')

p2 = ggplot(subset(df.mg,!is.na(TRS)),aes(x=as.factor(lin),y=scale(TRS))) + 
  geom_jitter(width=0.1,aes(col=scale(TRS)),alpha=0.6) + 
  geom_boxplot(outlier.shape=NA,alpha=0.8,fill='gray94',lwd=1,width=0.2,col="#FB8072") + 
  ggpubr::theme_pubr() +
  labs(x="HSC Branch",y=paste0("SCAVENGE TRS (",traitName,")")) +
  guides(col="none") +
  scale_color_gradient(high=viridis[length(viridis)],low=viridis[3]) + labs(col="") +
  scale_x_discrete(labels=paste0(tmp[,1],"\n (",round(tmp[,2]*100,1),"%)")); p2

