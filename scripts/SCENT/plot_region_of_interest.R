library(Seurat)
library(Signac)
library(ggplot2)
library(data.table)

dfseurat = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.HSC_only.rds")

lineage = fread("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Lineage_Probabilities_per_Cell.csv",data.table = F,stringsAsFactors = F)
mapping = fread("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/dataset_metadata_to_cell.csv",data.table = F,stringsAsFactors = F,header=T)
which(lineage$V1!=mapping$V1)
lineage = cbind(mapping,lineage[,-1])
lineage$cell = unlist(lapply(strsplit(lineage$V1,"_"),function(x) x[1]))
lineage$lin = apply(lineage[,paste0("Lineage_",c(1:3),"_prob")],1,which.max)
lineage$lin_prob_max = apply(lineage[,paste0("Lineage_",c(1:3),"_prob")],1,max)
lineage$lin[lineage$lin_prob_max < 0.8] = "NA"

dfseurat@meta.data$i = 1:nrow(dfseurat@meta.data)
df.mg = merge(dfseurat@meta.data,lineage,by=c("cell","dataset"),all=T)
df.mg = df.mg[order(df.mg$i),]

df.mg$lin2 = paste0(df.mg$lin,"-",df.mg$disease)

dfseurat@meta.data = df.mg

# gene_of_interest = "TFR2"; peak_of_interest = "chr7-100639909-100642992"
gene_of_interest = "TSPAN32"; peak_of_interest = "chr11-2300251-2303115"
Idents(dfseurat) <- "disease"
levels(dfseurat) = c("H","T21")
DefaultAssay(dfseurat) = "ATAC"
col_to_use = c("#FB8072" ,"#CBD5E8")
g = CoveragePlot(
  object = dfseurat,
  region = peak_of_interest,
  features = gene_of_interest,
  extend.upstream = 100,
  extend.downstream = 100
) & scale_fill_manual(values=col_to_use)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/plots/peak_coverage_plot.",gene_of_interest,".pdf")
pdf(f.plot,height=3.5,width=10)
print(g)
dev.off()
