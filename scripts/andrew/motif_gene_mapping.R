library(data.table)
library("JASPAR2020")
library("TFBSTools")
library(parallel)

res.all = fread("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_h_v2/h.ChromVAR.txt",data.table = F,stringsAsFactors = F)
# res.all = fread("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/ds.ChromVAR.txt",data.table = F,stringsAsFactors = F)
colnames(res.all)[1] = "gene"

print("getMatrixByID...")
genelst <- getMatrixByID(JASPAR2020, ID = unique(res.all$gene))

print("motif_name...")
motif_name <- as.character(
  unlist(
    mclapply(
      genelst,
      function(x) ID(x),
      mc.cores = 4)
  )
)

print("gene_name...")
gene_name <- as.character(
  unlist(
    mclapply(
      genelst,
      function(x) name(x),
      mc.cores = 4)
  )
)

print("motif-to-gene linking...")
res.df = data.frame(motif_name,gene_name)
f.out = "/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/motif_gene.ChromVAR.txt"
fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


