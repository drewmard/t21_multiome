for (k in 1:2) {
  
  
  if (k==1) {
    id="D21"
    df.sub=subset(data_snv,genotype==id)
  } else if (k==2) {
    id="T21"
    df.sub=subset(data_snv,genotype==id)
  }
  
  ind = !duplicated(df.sub$snv); sum(!ind); df.sub=df.sub[ind,]
  print(paste0(k," - ",id))
  
  snv.df=df.sub[,1:3]
  colnames(snv.df)[1:3] = c("chrom","start","end")
  snv.df$posid = paste0("snv",1:nrow(snv.df))
  snv.df$genotype = 1
  snv.df$strand = "+"
  
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/",id,".bed")
  fwrite(snv.df,f.out,col.names = F,row.names = F,quote = F,na = "NA",sep = '\t')
}

t21=/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/T21.bed
d21=/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/D21.bed
mkdir -p /oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/out
findMotifsGenome.pl $t21 hg38 /oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/out -size given -bg $d21

annotatePeaks.pl $t21 hg38

/oak/stanford/groups/smontgom/amarder/data/homer/homer.KnownMotifs.hg38.191020.bed.gz


