library(data.table)

celltype_to_use="HSCs_H"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
res.df.h = fread(output_file,data.table = F,stringsAsFactors = F)

celltype_to_use="HSCs_T21"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",celltype_to_use,".txt")
res.df.t21 = fread(output_file,data.table = F,stringsAsFactors = F)

res.df.mg = merge(res.df.h,res.df.t21,all=TRUE,by=c("gene","peak"))

cor.test(res.df.mg$beta.x,res.df.mg$beta.y,use='na.or.complete')
cor.test(res.df.mg$beta.x,res.df.mg$beta.y,use='na.or.complete')$p.value
cor.test(subset(res.df.mg,fdr.x<0.1 & fdr.y<0.1)$beta.x,subset(res.df.mg,fdr.x<0.1 & fdr.y < 0.1)$beta.y,use='na.or.complete')
cor.test(subset(res.df.mg,fdr.x<0.1)$beta.x,subset(res.df.mg,fdr.x<0.1)$beta.y,use='na.or.complete')

fisher.test(res.df.mg$fdr.x < 0.1, res.df.mg$fdr.y < 0.1)
res.df.mg[order(res.df.mg$fdr.y)[1:5],]

ggplot(subset(res.df.mg,!is.na(beta.x) & !is.na(beta.y)),aes(x=beta.x,y=beta.y)) + 
  geom_point() + 
  geom_abline(slope=1,intercept = 0,lty='dashed',col='red') +
  ggpubr::theme_pubr() +
  # theme_bw() +
  # theme(panel.grid = element_blank()) +
  # geom_smooth(se=F,col='red') +
  # geom_smooth(method='lm',se=F,col='red') +
  labs(x=bquote("Peak-gene"~beta~"(HSCs in control samples)"),
       y=bquote("Peak-gene"~beta~"(HSCs in T21 samples)")) 

subset(res.df.mg,fdr.x<0.1 & fdr.y<0.1)

ggplot(subset(res.df.mg,!is.na(beta.x) & !is.na(beta.y) & fdr.x<0.1 & fdr.y<0.1),aes(x=beta.x,y=beta.y)) + 
  geom_point() + 
  geom_abline(slope=1,intercept = 0,lty='dashed',col='red') +
  ggpubr::theme_pubr() +
  # theme_bw() +
  # theme(panel.grid = element_blank()) +
  # geom_smooth(se=F,col='red') +
  # geom_smooth(method='lm',se=F,col='red') +
  labs(x=bquote("Peak-gene"~beta~"(HSCs in control samples)"),
       y=bquote("Peak-gene"~beta~"(HSCs in T21 samples)")) 


ggplot(subset(res.df.mg,fdr.x < 0.1 & fdr.y < 0.1),aes(x=beta.x,y=beta.y)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') + theme_bw()
ggplot(subset(res.df.mg,fdr.y < 0.1),aes(x=beta.x,y=beta.y)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') + theme_bw()

# trait = fread("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/finemap/finemappedtraits_hg19/hg38/rbc.PP001.hg38.bed",data.table = F,stringsAsFactors = F)
trait_fm = fread("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/rbc_peaks_overlap.bed",data.table = F,stringsAsFactors = F)
trait_fm.sub = trait_fm[,c("V1","V3","V4","V5","V9")]
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")

# df.sub = subset(res.df.mg,peak %in% subset(trait_fm,V5 > 0.1)$V9 & 
#          (fdr.x < 0.1 | fdr.y < 0.1)
#        )
# ggplot(df.sub,aes(x=beta.x,y=beta.y)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') + theme_bw()
trait_fm.sub = subset(trait_fm,V5 > 0.2)[,c("V1","V3","V4","V5","V9")]
dim(trait_fm.sub)
nrow(subset(trait_fm,V5 > 0.9)[,c("V1","V3","V4","V5","V9")])
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")
merge(subset(res.df.mg,)

df.sub = merge(subset(res.df.mg,!is.na(beta.x)),subset(trait_fm.sub,pip>0.2),by="peak")
nrow(df.sub)
df.sub = merge(subset(res.df.mg,!is.na(beta.x)),subset(trait_fm.sub,pip>0.9),by="peak")
nrow(df.sub)
df.sub = merge(subset(res.df.mg,!is.na(beta.y)),subset(trait_fm.sub,pip>0.2),by="peak")
nrow(df.sub)
df.sub = merge(subset(res.df.mg,!is.na(beta.y)),subset(trait_fm.sub,pip>0.9),by="peak")
nrow(df.sub)

df.sub = merge(subset(res.df.mg,!is.na(beta.x) & fdr.x < 0.1),subset(trait_fm.sub,pip>0.2),by="peak")
nrow(df.sub)
df.sub = merge(subset(res.df.mg,!is.na(beta.x) & fdr.x < 0.1),subset(trait_fm.sub,pip>0.9),by="peak")
nrow(df.sub)
df.sub = merge(subset(res.df.mg,!is.na(beta.y) & fdr.y < 0.1),subset(trait_fm.sub,pip>0.2),by="peak")
nrow(df.sub)
df.sub = merge(subset(res.df.mg,!is.na(beta.y) & fdr.y < 0.1),subset(trait_fm.sub,pip>0.9),by="peak")
nrow(df.sub)

ggplot(df.sub,aes(x=beta.x,y=beta.y)) + geom_point() + geom_abline(slope=1,intercept = 0,lty='dashed',col='red') + theme_bw()
df.sub$amplify = log2(df.sub$beta.y/df.sub$beta.x)

df.sub = merge(subset(res.df.mg,!is.na(beta.y) & fdr.y < 0.1),subset(trait_fm.sub,pip>0.2),by="peak")
df.sub = df.sub[,c("chr","snp_pos","rsid","pip","peak","gene","beta.x","p.x","fdr.x","beta.y","p.y","fdr.y")]
colnames(df.sub)[7:12] = c("beta_h","p_h","fdr_h","beta_t21","p_t21","fdr_t21")
write.csv(df.sub,"~/Downloads/hsc_scent.csv")

da_peaks = fread("~/Downloads/da_peaks.hsc.txt",data.table = F,stringsAsFactors = F)
de_genes = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
df.sub2 = merge(df.sub,da_peaks,by.x='peak',by.y='id')
df.sub2 = merge(df.sub2,de_genes,by.x='gene',by.y='names')

head(df.sub2)

df.sub = merge(subset(res.df.mg,!is.na(beta.y) & fdr.y < 0.1),subset(trait_fm.sub,pip>0.2),by="peak")
nrow(df.sub)
ggplot(df.sub,aes(x=beta.x,y=beta.y)) + 
  geom_point() + 
  geom_abline(slope=1,intercept = 0,lty='dashed',col='red') +
  ggpubr::theme_pubr() +
  labs(x=bquote("Peak-gene"~beta~"(HSCs in control samples)"),
       y=bquote("Peak-gene"~beta~"(HSCs in T21 samples)")) +
  geom_vline(xintercept = 0,lty='dashed') +
  geom_hline(yintercept = 0,lty='dashed')
cor.test(df.sub$beta.x,df.sub$beta.y)
