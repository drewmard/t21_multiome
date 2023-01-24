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
gene_of_interest = "ABCA7"; peak_of_interest = "chr19-1063288-1063650"
# gene_of_interest = "KLF6"; peak_of_interest = "chr10-3774036-3776511"
# gene_of_interest = "HBS1L"; peak_of_interest = "chr6-135096824-135097623"
gene_of_interest = "FNTB"; peak_of_interest = "chr14-65042220-65043484"
gene_of_interest = "TMC6"; peak_of_interest = "chr17-78128031-78129288"
gene_of_interest = "TSPAN32"; peak_of_interest = "chr11-2300251-2303115"
# Idents(dfseurat) <- "lin2"
# levels(dfseurat) = c("1-H","1-T21","2-H","2-T21","3-H","3-T21")
Idents(dfseurat) <- "disease"
levels(dfseurat) = c("H","T21")
DefaultAssay(dfseurat) = "ATAC"
col_to_use = c("#FB8072" ,"#CBD5E8")
# col_to_use = c("#FB8072" ,"#CBD5E8","violetred3")
# col_to_use = c("#FB8072","#FB8072" ,"#CBD5E8","#CBD5E8","violetred3","violetred3")
g = CoveragePlot(
  object = dfseurat,
  region = peak_of_interest,
  features = gene_of_interest,
  extend.upstream = 100,
  extend.downstream = 100
) & scale_fill_manual(values=col_to_use)
f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/plots/peak_coverage_plot.",gene_of_interest,".pdf")
# pdf(f.plot,height=8,width=10)
pdf(f.plot,height=3.5,width=10)
print(g)
dev.off()
