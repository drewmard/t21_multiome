# Run on cluster!
# Needed for input into Fig K.

dfseurat.sub.rna = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA.HSC_only.rds")
# dfseurat.sub.atac = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

traitName = "rbc"
trait_mat.save = fread(paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/old/scavenge/",traitName,".txt"),data.table = F,stringsAsFactors = F)
hsc.sub = subset(trait_mat.save,subclust_v6=="HSCs" & TRS > quantile(trait_mat.save$TRS,probs=0.9))
dfseurat.sub.rna$rbc = dfseurat.sub.rna$cell %in% hsc.sub$cell

dfseurat.sub.rna
gene_of_interest_lst = c("GATA1","ITGA2B","MKI67","PROM1","SPINK2")
for (gene_of_interest in gene_of_interest_lst) {
  print(gene_of_interest)
  tmp = data.frame(trs=dfseurat.sub.rna$rbc,expr=as.numeric(dfseurat.sub.rna@assays$RNA@data[gene_of_interest,]))
  print(aggregate(.~trs,tmp,median))
  print(aggregate(.~trs,tmp,mean))
}


expr.df = apply(t(dfseurat.sub.rna@assays$RNA@data[gene_of_interest_lst,]),2,as.numeric)
# expr.df=data.frame(cell=dfseurat.sub.rna$cell,trs_ex=dfseurat.sub.rna$rbc,expr.df)
expr.df=data.frame(cell=dfseurat.sub.rna$cell,disease=dfseurat.sub.rna$disease,trs_ex=dfseurat.sub.rna$rbc,expr.df)
fwrite(expr.df,"/home/amarder/tmp/expr.df.txt",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
