library(data.table)
library(Matrix.utils)
library(Seurat)
library(DESeq2)
library(variancePartition)
library('edgeR')

print("Reading metadata2...")
disease_status="DownSyndrome"

cell_type="Cycling HSCs/MPPs";
disease_status="DownSyndrome"
sampletype="Liver"

pseudobulk_data <- function(cell_type,disease_status,sampletype) {
  
  print(paste(cell_type,disease_status,sampletype))
  print(paste0("Metadata..."))
  f = paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/cellComp/10X_",disease_status,"_",sampletype,".cellComp.csv")
  meta2.full<-fread(f,data.table = F,stringsAsFactors = F)
  colnames(meta2.full)[6] <- "leiden_names"
  cells2 <- unique(meta2.full[,6])
  meta2 <- unique(meta2.full[,c("patient","sample","sorting")])
  meta2$environment <- disease_status
  meta2$sorting[!(meta2$sorting %in% c("CD235a-","CD45+"))] <- "Other"
  x <- meta2
  rownames(x) <- x[,"sample"]
  rm(meta2)
  
  print(paste0("Cell data..."))
  
  cell_type_filename = gsub("/","_",cell_type)
  f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
  df <- readRDS(file = f.out)
  
  print(paste0("Pseudobulking..."))
  
  df.aggre <- aggregate.Matrix(
    t(
      GetAssayData(object = df, slot = "counts", assay="RNA")
    ),
    groupings=df@meta.data[,"sample"],fun="sum")
  df.aggre <- t(df.aggre)
  df.aggre <- as.data.frame(df.aggre)
  
  print(paste0("Processing..."))
  
  metadata_to_use <- x[rownames(x) %in% colnames(df.aggre),]
  df.aggre <- as.matrix(df.aggre[,match(rownames(metadata_to_use),colnames(df.aggre))])
  
  tab <- table(meta2.full[meta2.full[,"leiden_names"]==cell_type,'sample'])
  samples_to_keep <- names(tab)[tab >= 10]
  samples_to_keep <- samples_to_keep[samples_to_keep %in% rownames(metadata_to_use)]
  df.aggre <- df.aggre[,samples_to_keep]
  metadata_to_use <- metadata_to_use[samples_to_keep,]
  df.aggre1=df.aggre
  metadata_to_use1 <- metadata_to_use
  metadata_to_use1$leiden_names=cell_type
  if (cell_type=="Cycling HSCs/MPPs") {
    colnames(df.aggre1) <- paste0(colnames(df.aggre1),".cyc")
    rownames(metadata_to_use1) <- paste0(rownames(metadata_to_use1),".cyc")
  }
  
  df.aggre1 = as.data.frame(df.aggre1)
  
  return(list(df.aggre1,metadata_to_use1))
}

cell_type="Cycling HSCs/MPPs";
disease_status="DownSyndrome"
sampletype="Liver"
pb_lst = pseudobulk_data(cell_type=cell_type,disease_status=disease_status,sampletype=sampletype)
# dim(pb_lst[[2]])

cell_type_filename = gsub("/","_",cell_type)
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_expr.txt")
fwrite(pb_lst[[1]],f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)

f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_meta.txt")
fwrite(pb_lst[[2]],f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)


cell_type="HSCs/MPPs";
disease_status="DownSyndrome"
sampletype="Liver"
pb_lst = pseudobulk_data(cell_type=cell_type,disease_status=disease_status,sampletype=sampletype)
# dim(pb_lst[[2]])

cell_type_filename = gsub("/","_",cell_type)
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_expr.txt")
fwrite(pb_lst[[1]],f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)

f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_meta.txt")
fwrite(pb_lst[[2]],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

cell_type="HSCs/MPPs";
disease_status="Healthy"
sampletype="Liver"
pb_lst = pseudobulk_data(cell_type=cell_type,disease_status=disease_status,sampletype=sampletype)
dim(pb_lst[[2]])

cell_type_filename = gsub("/","_",cell_type)
f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_expr.txt")
fwrite(pb_lst[[1]],f.out,quote = F,na = "NA",sep = '\t',row.names = T,col.names = T)

f.out=paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_meta.txt")
fwrite(pb_lst[[2]],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)


########################################

library(data.table)
library(edgeR)

read_data_calc_cpm <- function(cell_type,disease_status,sampletype,calc_cpm=FALSE) {
  cell_type_filename = gsub("/","_",cell_type)
  f=paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_expr.txt")
  df1=fread(f,data.table = F,stringsAsFactors = F)
  f=paste0("/Users/andrewmarderstein/Documents/Research/t21-proj/out/full/pb/",sampletype,".",disease_status,".",cell_type_filename,".pb_meta.txt")
  meta1=fread(f,data.table = F,stringsAsFactors = F)
  rownames(df1) = df1[,1]
  df1 = df1[,-1]
  if (calc_cpm) {df1 = cpm(df1)}
  return(list(df1,meta1))
}

df = list()
cell_type="Cycling HSCs/MPPs";
disease_status="DownSyndrome"
sampletype="Liver"
df[[1]] = read_data_calc_cpm(cell_type,disease_status,sampletype)
cell_type="HSCs/MPPs";
disease_status="DownSyndrome"
sampletype="Liver"
df[[2]] = read_data_calc_cpm(cell_type,disease_status,sampletype)
cell_type="HSCs/MPPs";
disease_status="Healthy"
sampletype="Liver"
df[[3]] = read_data_calc_cpm(cell_type,disease_status,sampletype,calc_cpm=FALSE)

df.aggre = merge(df[[1]][[1]],df[[2]][[1]],by='row.names')
df.aggre = merge(df.aggre,df3,by.x="Row.names",by.y='row.names')
rownames(df.aggre) = df.aggre[,1]
df.aggre = df.aggre[,-1]
df.aggre.cpm = cpm(df.aggre)

meta = rbind(df[[1]][[2]][,-1],
      df[[2]][[2]],
      df[[3]][[2]])
meta$group = NA
meta$group[meta$leiden_names=="Cycling HSCs/MPPs"] <- "Ts21 cycling HSCs"
meta$group[meta$leiden_names=="HSCs/MPPs" & meta$environment=="DownSyndrome"] <- "Ts21 HSCs"
meta$group[meta$leiden_names=="HSCs/MPPs" & meta$environment=="Healthy"] <- "Healthy HSCs"

geneName="TFR2"
out = data.frame(group=meta$group,expr=as.numeric(df.aggre.cpm[geneName,]))
out$group = factor(out$group, rev(c("Healthy HSCs", "Ts21 HSCs","Ts21 cycling HSCs")))
col_to_use = c("darkred" ,"#CBD5E8","#FB8072")

# aggregate(expr~group,out,mean)
g=ggplot(out,aes(x=group,y=expr,fill=group)) +
  theme_bw() +
  # geom_violin() +
  geom_boxplot(outlier.shape=NA,width=0.2) +
  geom_point(alpha=0.6) +
  scale_fill_manual(values=col_to_use) +
  # scale_fill_brewer(palette = "Set1") +
  theme(panel.grid = element_blank()) + 
  labs(x="Group",y=bquote(italic(.(geneName))~"expression in large scRNA-seq (cpm)")) + coord_flip() +
  guides(fill="none") +
  scale_y_continuous(trans = "log10") +
  ggpubr::theme_pubr(); g

f.plot = paste0("~/Documents/Research/t21_multiome/output/scent/plots/HSC_pb.",geneName,".cyc.pdf")
pdf(f.plot,height=3.1,width=5.3)
print(g)
dev.off()


