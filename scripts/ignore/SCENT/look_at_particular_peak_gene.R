
celltype_to_use="HSCs_t21"
gene="PEX14"
this_peak="chr1-10474439-10475771"

# meta data:
meta = data.frame(
  cell=rownames(dfseurat@meta.data),
  percent_mito = dfseurat@meta.data$percent.mt,
  nUMI=colSums(dfseurat@assays$RNA@counts),
  sample=dfseurat@meta.data$dataset,
  celltype0 = gsub(" ","_",dfseurat@meta.data$subclust_v6),
  disease0 = dfseurat@meta.data$disease,
  celltype=dfseurat@meta.data$celltype
)
rownames(meta) <- NULL

atac_df<-data.frame(cell=colnames(dfseurat@assays$ATAC@data),
                        atac_data=as.numeric(dfseurat@assays$ATAC@data[this_peak,]),
                        atac_count=as.numeric(dfseurat@assays$ATAC@counts[this_peak,]))
expr_df <- data.frame(cell=colnames(dfseurat@assays$RNA@data),
                 expr_data=as.numeric(dfseurat@assays$RNA@data[gene,]),
                 expr_count=as.numeric(dfseurat@assays$RNA@counts[gene,]))
df<-merge(atac_df,expr_df,by="cell")
df<-merge(df,meta,by="cell")

# Subset cells to test:
celltype_to_use="HSCs_T21"
df2 <- df[df$celltype==celltype_to_use,]
celltype_to_use="HSCs"
df2 <- df[df$celltype0==celltype_to_use,]
aggregate(expr_data ~ atac_data + percent_mito + nUMI,data=df2,mean)
aggregate(expr_data ~ disease0 + (atac_data > 0),data=df2,mean)

summary(lm(expr_data ~ atac_data + percent_mito + nUMI + sample,data=df))
summary(lm(expr_data ~ atac_data + percent_mito + nUMI + sample,data=subset(df,celltype=="HSCs_H")))
summary(lm(expr_data ~ atac_data + percent_mito + nUMI + sample,data=subset(df,celltype=="HSCs_T21")))
summary(lm(expr_data ~ (atac_data > 0) + percent_mito + nUMI + sample,data=subset(df,celltype=="HSCs_T21")))
summary(lm(expr_count ~ (atac_data > 0) + percent_mito + nUMI + sample,data=subset(df,celltype=="HSCs_T21")))
summary(lm(expr_count ~ (atac_data > 0) + percent_mito + nUMI,data=subset(df,celltype=="HSCs_T21")))
summary(lm(expr_data ~ atac_data + percent_mito + nUMI,data=subset(df,celltype=="HSCs_H")))
summary(lm(expr_data ~ atac_data*disease0 + percent_mito + nUMI + sample,data=subset(df,celltype%in%c("HSCs_H","HSCs_T21"))))

summary(lm(expr_data ~  (atac_data > 0),subset(df,disease0=="T21" & celltype0=="HSCs")))
summary(lm(expr_data ~ (atac_data > 0),data=subset(df,celltype=="HSCs_T21")))
summary(lm(expr_data ~ (atac_data > 0) + percent_mito + sample,data=subset(df,celltype=="HSCs_T21")))
summary(lm(expr_data ~ (atac_data > 0) + percent_mito + sample,data=subset(df,celltype=="HSCs_H")))
summary(lm(expr_data ~ disease0,data=subset(df,celltype0=="HSCs")))


summary(lm(expr_data ~ (atac_data > 0) + percent_mito + nUMI + sample,data=subset(df,celltype=="HSCs_T21")))


# Binarize peaks:
df2[df2$atac>0,]$atac<-1
