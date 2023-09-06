# snv
k=2
if (k==1) {
  id="D21"
  df.sub=subset(data_snv,genotype==id)
} else if (k==2) {
  id="T21"
  df.sub=subset(data_snv,genotype==id)
}

ind = !duplicated(df.sub$snv); sum(!ind); df.sub=df.sub[ind,]
print(paste0(k," - ",id))

snv.df=df.sub
colnames(snv.df)[1:3] = c("chrom","start","end")

# snvs in genes
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror = "useast")
annot <- getBM(attributes = c('hgnc_symbol','chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol',
               values = atac$names,
               mart = ensembl)
annot = subset(annot,chromosome_name %in% c(1:22,"X","Y"))
annot = annot[!duplicated(annot$hgnc_symbol),]
colnames(annot)[2:4] = c("chrom","start","end")
annot$chrom = paste0("chr",annot$chrom)

df.mg = as.data.frame(valr::bed_intersect(snv.df,annot))
df.mg = df.mg[,c("snv.x","hgnc_symbol.y")]
colnames(df.mg) = c("snv","hgnc_symbol")
snv.df = merge(snv.df,df.mg,by="snv",all.x=TRUE)

# snvs in exons
annot <- getBM(attributes = c('hgnc_symbol','ensembl_gene_id'),
               filters = 'hgnc_symbol', 
               values = atac$names, 
               mart = ensembl)
annot2 <- getBM(attributes = c('ensembl_gene_id','chromosome_name',
                               'exon_chrom_start',"exon_chrom_end"),
                filters = 'ensembl_gene_id', 
                values = annot$ensembl_gene_id, 
                mart = ensembl)

annot2 = subset(annot2,chromosome_name %in% c(1:22,"X","Y"))
annot = merge(annot,annot2,by="ensembl_gene_id")
colnames(annot)[3:5] = c("chrom","start","end")
annot$chrom = paste0("chr",annot$chrom)

snv.df.exon = as.data.frame(valr::bed_intersect(annot,snv.df))
snv.df.exon2 = snv.df.exon[,c("hgnc_symbol.x","chrom","start.y","end.y","fetus.y","genotype.y","snv.y")]
colnames(snv.df.exon2) = c("hgnc_symbol","chrom","start","end","fetus","genotype","snv")
ind = !duplicated(snv.df.exon2$snv); sum(!ind); snv.df.exon2=snv.df.exon2[ind,]
snv.df$exon = snv.df$snv %in% snv.df.exon2$snv
snv.df$gene = !is.na(snv.df$hgnc_symbol)
snv.df$intron = snv.df$gene & !snv.df$exon

# snv analysis:
snv.df.use= subset(snv.df,gene)
snv.df.use= subset(snv.df,intron)
atac.df=atac
colnames(atac.df)[8:10] = c("chrom","start","end")
atac.df$chrom = paste0("chr",atac.df$chrom)

df.mg = as.data.frame(valr::bed_intersect(atac.df,snv.df.use))
tab = table(df.mg$names.x)
tab = as.data.frame(tab)

atac.df$snv = atac.df$names %in% df.mg$names.x
atac.df = merge(atac.df,tab,by.x="names",by.y="Var1",all.x=TRUE)
atac.df$Freq[is.na(atac.df$Freq)] = 0
atac.df$span = atac.df$end - atac.df$start

bootstrap_snv_rate = function(j) {
  {i=sample(1:nrow(tmp),nrow(tmp),replace = T); tmp2 = tmp[i,]; val = sum(tmp2$Freq)/sum(tmp2$span); return(val)}
}

iter = 0; res = list()

iter=iter+1
tmp = subset(atac.df,adj.P.Val < 0.01 & logFC < -1)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:1000,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Cycling HSC downregulated (FDR < 0.01; LFC < -1)",est=val.open,l=ci[1],h=ci[2])

iter=iter+1
tmp = subset(atac.df,adj.P.Val < 0.01 & logFC > 1)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:1000,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Cycling HSC upregulated (FDR < 0.01; LFC > 1)",est=val.open,l=ci[1],h=ci[2])

iter=iter+1
tmp = subset(atac.df,adj.P.Val > 0.01 & logFC < 0.25 & logFC > -0.25)
# tmp = subset(atac.df,adj.P.Val > 0.01) #logFC < 1 & logFC > -1)
val.open = sum(tmp$Freq)/sum(tmp$span)
boot.open = unlist(parallel::mclapply(1:1000,bootstrap_snv_rate,mc.cores = 8))
ci = quantile(boot.open,probs=c(0.025,0.975))
res[[iter]] = data.frame(type="Background (FDR > 0.01; -0.25 < LFC < 0.25)",est=val.open,l=ci[1],h=ci[2])
# 
out.df = as.data.frame(do.call(rbind,res))
print(out.df)

