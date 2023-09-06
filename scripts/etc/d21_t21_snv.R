library(data.table)

# metadata
meta = fread("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/metadata.csv",data.table = F,stringsAsFactors = F)

# for loop around all files
headdir="/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38"
files = list.files(path=headdir,full.names = TRUE)
files_without_string = files[!grepl("unmapp",files)]
files_without_string = files_without_string[grep("bed",files_without_string)]
iter=0; df=list()
for (inFile in files_without_string) {
  
  iter = iter + 1
  
  # read in mutation data
  tmp=fread(inFile,data.table = F,stringsAsFactors = F)
  
  # add fetus and genotype info
  fetusName=sub("_complete.bed","",sub(paste0(headdir,"/"),"",inFile))
  genotypeInfo=meta[meta$name==fetusName,"genotype"]
  
  tmp$fetus = fetusName
  tmp$genotype = genotypeInfo
  
  # save data
  df[[iter]] = tmp
  
  # end for loop
}

# merge together
data_snv = as.data.frame(do.call(rbind,df))
data_snv$snv = paste(data_snv[,1],data_snv[,2],data_snv[,3],sep = '-')
fwrite(data_snv,"/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/data_snv",quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)

##############

disease_status="DownSyndrome"
pseudobulk_covariate="sample"
suffix=""
f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pseudobulks/",pseudobulk_covariate,"/DE/",disease_status,'_',"cyc_vs_hsc",".de",suffix,".txt")
atac=fread(f,data.table = F,stringsAsFactors = F)[,1:7]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror = "useast")
annot <- getBM(attributes = c('hgnc_symbol','chromosome_name',
                              'start_position', 'end_position'),
               filters = 'hgnc_symbol',
               values = atac$names,
               mart = ensembl)
annot = subset(annot,chromosome_name %in% c(1:22,"X","Y"))
annot = annot[!duplicated(annot$hgnc_symbol),]
atac = merge(atac,annot,by.x="names",by.y="hgnc_symbol")

for (k in 1:2) {

  if (k==1) {
    id="D21"
    df.sub=subset(data_snv,genotype==id)
  } else if (k==2) {
    id="T21"
    df.sub=subset(data_snv,genotype==id)
  } else if (k==3) {
    id="T21_MyeloidPreleuk"
    df.sub=subset(data_snv,genotype==id)
  } else if (k==4) {
    id="T21_and_T21_MyeloidPreleuk"
    df.sub=subset(data_snv,genotype%in%c("T21_MyeloidPreleuk","T21"))
  }
  ind = !duplicated(df.sub$snv); sum(!ind); df.sub=df.sub[ind,]
  
  
  print(paste0(k," - ",id))
  
  snv.df=df.sub
  colnames(snv.df)[1:3] = c("chrom","start","end")
  
  atac.df=atac
  colnames(atac.df)[8:10] = c("chrom","start","end")
  atac.df$chrom = paste0("chr",atac.df$chrom)

  df.mg = as.data.frame(valr::bed_intersect(atac.df,snv.df))
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
  tmp = subset(atac.df,adj.P.Val < 0.05 & logFC < 0)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:1000,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  res[[iter]] = data.frame(type="Cycling HSC downregulated (FDR < 0.01; LFC < -1)",est=val.open,l=ci[1],h=ci[2])
  
  iter=iter+1
  # tmp = subset(atac.df,adj.P.Val < 0.01 & logFC > 1)
  tmp = subset(atac.df,adj.P.Val < 0.05 & logFC > 0)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:1000,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  res[[iter]] = data.frame(type="Cycling HSC upregulated (FDR < 0.01; LFC > 1)",est=val.open,l=ci[1],h=ci[2])
  
  iter=iter+1
  tmp = subset(atac.df,adj.P.Val > 0.1 & logFC < 0.5 & logFC > -0.5)
  val.open = sum(tmp$Freq)/sum(tmp$span)
  boot.open = unlist(parallel::mclapply(1:1000,bootstrap_snv_rate,mc.cores = 8))
  ci = quantile(boot.open,probs=c(0.025,0.975))
  data.frame(type="Background (FDR > 0.01; -0.25 < LFC < 0.25)",est=val.open,l=ci[1],h=ci[2])
  res[[iter]] = data.frame(type="Background (FDR > 0.01; -0.25 < LFC < 0.25)",est=val.open,l=ci[1],h=ci[2])
  # 
  out.df = as.data.frame(do.call(rbind,res))
  print(out.df)
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/out/",id,".cycling.snv_enrich.txt")
  # fwrite(out.df,f.out,quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
  
}





##########################

library(data.table)
id="D21"
f=paste0("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/",id,".intron.cre.bed")
intron_cre=fread(f,data.table = F,stringsAsFactors = F)
f=paste0("/oak/stanford/groups/smontgom/amarder/data/hasaart_et_al_2020_scirep/hg38/homer/",id,".intron.bed")
intron=fread(f,data.table = F,stringsAsFactors = F)
dim(intron)
dim(intron_cre)

atac.df=atac
colnames(atac.df)[8:10] = c("chrom","start","end")
atac.df$chrom = paste0("chr",atac.df$chrom)

df.mg = as.data.frame(valr::bed_intersect(atac.df,snv.df))
tab = table(df.mg$names.x)
tab = as.data.frame(tab)

atac.df$snv = atac.df$names %in% df.mg$names.x
atac.df = merge(atac.df,tab,by.x="names",by.y="Var1",all.x=TRUE)
atac.df$Freq[is.na(atac.df$Freq)] = 0
atac.df$span = atac.df$end - atac.df$start

