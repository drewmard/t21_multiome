# snv
# snvs in genes
# module load R/4.1.2
library(biomaRt)
library(data.table)
library(valr,lib.loc="/home/amarder/bin")

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror = "useast")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

cre = fread("/oak/stanford/groups/smontgom/amarder/data/GRCh38-cCREs.bed",data.table=F,stringsAsFactors = F)
# cre = fread("/oak/stanford/groups/smontgom/amarder/data/ENCFF394DBM.bed",data.table=F,stringsAsFactors = F)
colnames(cre)[1:3] = c("chrom",'start','end')

for(k in 1:2) { 
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
  snv.df= subset(snv.df,intron)
  # snv.df.use= subset(snv.df,intron)
  
  snv.df=snv.df[,c("chrom","start","end")]
  # snv.df=df.sub[,1:3]
  colnames(snv.df)[1:3] = c("chrom","start","end")
  snv.df$posid = paste0("snv",1:nrow(snv.df))
  snv.df$genotype = 1
  snv.df$strand = "+"
  # snv.df$start=snv.df$start-10
  # snv.df$end=snv.df$end+10
  
  df.mg = as.data.frame(valr::bed_intersect(cre,snv.df))
  print(nrow(df.mg))
  print(nrow(snv.df))
  
  # print(nrow(snv.df))
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/",id,".intron.bed")
  fwrite(snv.df,f.out,col.names = F,row.names = F,quote = F,na = "NA",sep = '\t')
  
  f.out=paste0("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/",id,".intron.cre.bed")
  fwrite(df.mg,f.out,col.names = F,row.names = F,quote = F,na = "NA",sep = '\t')
}


atac.df=merge(atac,df.mg,by.x="names",by.y="hgnc_symbol.y")
library(data.table)
f.out="/oak/stanford/groups/smontgom/amarder/tmp/t21.intron.cre.genes.txt"
fwrite(atac.df[,c("names","logFC","P.Value","adj.P.Val")],f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
