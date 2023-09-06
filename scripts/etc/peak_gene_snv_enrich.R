library(data.table)
e2g = fread("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out_split/all/HSCs.DE.int.txt",data.table = F,stringsAsFactors = F)
y=strsplit(e2g$peak,"-")
e2g$chrom = paste0(unlist(lapply(y,function(x)x[[1]])))
e2g$start = as.numeric(paste0(unlist(lapply(y,function(x)x[[2]]))))-1
e2g$end = as.numeric(paste0(unlist(lapply(y,function(x)x[[3]]))))

snv.df.keep=list()
for (k in 1:2) { 
  print(k)
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
  snv.df$genotype = id

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
  
  snv.df.keep[[k]] = snv.df
}

snv.df.all = do.call(rbind,snv.df.keep)

df.mg = as.data.frame(valr::bed_intersect(e2g,snv.df.all))
df.sub = subset(df.mg,intron.y)
# df.sub = df.mg
table(df.sub[,c("genotype.y")],df.sub[,"P.Value_smRNA.x"] < 0.05)
table(df.sub[,c("genotype.y")],df.sub[,"P.Value_smATAC.x"] < 0.05)
table(df.sub[,c("genotype.y")],df.sub[,"boot_basic_p_int.x"] < 0.05)
table(df.sub[,c("genotype.y")],df.sub[,"boot_basic_p_t21.x"] < 0.05)


small_ATAC_pb = fread("~/Downloads/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)
# repeat analysis here?

matrix()


