f.out = paste0("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/Multiome.RNA_ATAC.h_and_t21.ChromVAR.rds")
dfcombined=readRDS(f.out)
motif_gene = fread("/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/DS_Multiome_ds_v2/motif_gene.ChromVAR.txt",data.table = F,stringsAsFactors = F)
aggregate(as.numeric(dfcombined[['chromvar']]['MA0035.4',]),
          by=list(dataset=dfcombined@meta.data$dataNum),
          median
)
s=(unlist(lapply(strsplit(colnames(chromvar),"_"),function(x)x[[1]])))
table(s %in% dfseurat@meta.data$cell)
sum(colnames(chromvar) %in% colnames(dfseurat))
sum(colnames(chromvar) %in% colnames(dfseurat))

dfseurat.sub.atac = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

chromvar = dfcombined@assays$chromvar@data[,match(x2,x1)]
colnames(chromvar) = rownames(dfseurat.sub.atac@meta.data)
dfseurat.sub.atac[["chromvar"]] = CreateAssayObject(data = chromvar)

f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.atac_umap.GATA1.pdf")
pdf(f.plot,width=6,height=4)
DefaultAssay(dfseurat.sub.atac) <- "chromvar"
print(
  FeaturePlot(dfseurat.sub.atac,features=c("MA0035.4")) + scale_colour_gradient2(midpoint = 0,mid='grey')
)
# DefaultAssay(dfseurat.sub.rna) <- "ATAC"
dev.off()
aggregate(as.numeric(dfseurat.sub.atac[['chromvar']]['MA0035.4',]),
          by=list(disease=dfseurat.sub.atac@meta.data$disease),
          median
)
tmp = data.frame(gata1=as.numeric(dfseurat.sub.atac[['chromvar']]['MA0035.4',]),
           disease=dfseurat.sub.atac@meta.data$disease,
           kmeans=dfseurat.sub.atac@meta.data$)
library(ggplot2)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/HSCs.violin.GATA1.pdf")
pdf(f.plot,width=6,height=4)
print(ggplot(tmp,aes(x=disease,y=gata1)) + geom_violin())
dev.off()

which(dfcombined$cell_dataset %in% dfseurat.sub.atac$cell_dataset)

x1=dfcombined$cell_dataset
x2=dfseurat.sub.atac$cell_dataset
which(x1[match(x2,x1)]!=x2)
which(x1[match(x2,x1)]==x2)






