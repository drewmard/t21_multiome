###############################################################

# Load:
library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(gchromVAR)
library(SCAVENGE)
library(harmony)
library(dplyr)
library(uwot)
library(parallel)
library(ggplot2)
library(data.table)

###############################################################

# Input version 1:
mutualknn30 = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/shared_knn.rds")
dfseurat = readRDS("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/data/Multiome.RNA_ATAC.all.rds")

# Input version 2:
# input cell data needs to have:
# - scATAC assay called "ATAC"
# - harmony embeddings for KNN
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
dfseurat = readRDS(input_file)
mutualknn30 <- getmutualknn(dfseurat@reductions$harmony@cell.embeddings, 30)

# Set:
numThreads = detectCores()/2
viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF",
            "#5DC863FF", "#AADC32FF", "#FDE725FF")

####################################################

# transform seurat to summarizedexperiment
SE_data <- SummarizedExperiment(assays = list(counts = dfseurat@assays$ATAC@counts),
                                rowData = dfseurat@assays$ATAC@ranges, 
                                colData = DataFrame(names = colnames(dfseurat@assays$ATAC),dataset=dfseurat@meta.data$dataset))

SE_data <- addGCBias(SE_data, genome = BSgenome.Hsapiens.UCSC.hg38)

# preprocess gchromvar
SE_data_bg <- getBackgroundPeaks(SE_data, niterations=200)

# Pull out hg38 fine-mapped bed files:
fileDir="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/finemap/finemappedtraits_hg19/hg38"
flst = list.files(fileDir)
flst = list.files(fileDir,pattern = "hg38",include.dirs = FALSE)

# Cycling through every trait:
trait_mat.save.all = list()
for (k in 1:length(flst)) {
  
  # Pull trait info:
  f = flst[k]
  trait_file = paste0(fileDir,"/",f)
  traitName = gsub("\\..*","",f)
  print(paste0(k,"/",length(flst),": ",traitName))
  
  # Read in SNP data and perform gchromvar Z-scores:
  trait_import <- importBedScore(rowRanges(SE_data), trait_file, colidx=5)
  SE_data_DEV <- computeWeightedDeviations(SE_data, trait_import, background_peaks =
                                             SE_data_bg)
  z_score_mat <- data.frame(colData(SE_data), z_score=t(assays(SE_data_DEV)[["z"]]) %>% c)
  
  # Seed cells:
  seed_idx <- seedindex(z_score_mat$z_score, 0.05)
  scale_factor <- cal_scalefactor(z_score=z_score_mat$z_score, 0.01)
  
  # Network propagation from the seed cells
  np_score <- randomWalk_sparse(intM=mutualknn30, rownames(mutualknn30)[seed_idx], gamma=0.05)
  
  # remove cells not connected to other cells (singleton cells)
  omit_idx <- np_score==0
  mutualknn30.2 <- mutualknn30[!omit_idx, !omit_idx]
  np_score <- np_score[!omit_idx]
  
  # Calculate TRS:
  TRS <- np_score %>% capOutlierQuantile(., 0.95) %>% max_min_scale
  TRS <- TRS * scale_factor
  
  # Output matrix:
  trait_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx[!omit_idx], np_score, TRS)
  trait_mat = merge(data.frame(names=rownames(dfseurat@meta.data),
                               cell=dfseurat@meta.data$cell,
                               disease=dfseurat@meta.data$disease,
                               subclust_v6=dfseurat@meta.data$subclust_v6,
                               UMAP_1=dfseurat@reductions$umap@cell.embeddings[,1],
                               UMAP_2=dfseurat@reductions$umap@cell.embeddings[,2]),
                    trait_mat,
                    by='names')
  
  nperm = 1000
  trait_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30.2, seed_idx=trait_mat$seed_idx,
                                    topseed_npscore=trait_mat$np_score, permutation_times=nperm,
                                    true_cell_significance=0.05, rda_output=F, mycores=numThreads, rw_gamma=0.05)
  colnames(trait_permu)[2] <- "enriched"
  trait_mat <- data.frame(trait_mat, trait_permu)
  # trait_mat[order(trait_mat$TRS,decreasing = T)[1:30],]
  
  # Save output matrix:
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge/",traitName,".txt")
  trait_mat.save = trait_mat[,c("names","cell","dataset","disease","subclust_v6","UMAP_1","UMAP_2","TRS","enriched")]
  fwrite(trait_mat.save,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
  # TRS U-MAP:
  p2 <- ggplot(data=trait_mat.save, aes(UMAP_1, UMAP_2, color=TRS)) + 
    geom_point(size=1, na.rm = TRUE, alpha = 0.6) +
    scale_color_gradientn(colors = viridis) + 
    scale_alpha()+
    theme_bw() + 
    theme(legend.title = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
  f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_plots/",traitName,".pdf")
  pdf(f.plot,width=7,height=7)
  print(p2)
  dev.off()
  
  # TRS boxplot:
  g=ggplot(trait_mat.save,aes(x=subclust_v6,y=TRS,fill=disease)) +
    theme_bw() +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1") +
    theme(panel.grid = element_blank(),axis.text.x = element_text(angle=45,vjust=1,hjust=1)) +
    labs(x="Cell type",y=paste0("TRS (",traitName,")"))
  f.plot = paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scavenge_boxplots/boxplot.",traitName,".pdf")
  pdf(f.plot,width = 10,height=6)
  print(g)
  dev.off()
  
}






