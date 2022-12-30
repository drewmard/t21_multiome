library(data.table)

# pb_da_peaks = fread("~/Downloads/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)

pct.rna.keep = 0.05
pct.atac.keep = 0.05
dfseurat$kmeans_RNA = tmp$kmeans_RNA
atac.all <- dfseurat@assays$ATAC@counts
atac.all@x[atac.all@x > 1] <- 1

peak.data.lst = list()
for (i in 1:3) { 
  print(i)
  atac.all.sub = atac.all[,dfseurat$kmeans_RNA==i]
  pct.accessible = rowSums(atac.all.sub)/ncol(atac.all.sub)
  peak.data = data.frame(peakid=rownames(atac.all),pct=pct.accessible)
  peak.data.lst[[i]] = peak.data
}

peak.data.df = as.data.frame(do.call(cbind,peak.data.lst))[,c(1,2,4,6)]
rownames(peak.data.df) = NULL
colnames(peak.data.df)[2:4] = paste0("acc",1:3)

fwrite(peak.data.df,"/home/amarder/tmp/peak_df.txt",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

trait_peaks_all = fread("/home/amarder/tmp/rbc_peaks_overlap.bed",data.table = F,stringsAsFactors = F)
PIP_THRES=0.2; cat(paste0("\n\nThere are ",nrow(subset(trait_peaks_all,V5 > PIP_THRES))," fine-mapped SNPs (PIP>",PIP_THRES,") in ",length(snp_peaks <- unique(subset(trait_peaks_all,V5 > PIP_THRES)$V9))," peaks.\n\n"))
peak.df = subset(peak.data.df,peakid %in% snp_peaks)
fwrite(peak.df,paste0("/home/amarder/tmp/rbc.fm_peak.pip_",PIP_THRES,".txt"),quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)
PIP_THRES=0.9; cat(paste0("\n\nThere are ",nrow(subset(trait_peaks_all,V5 > PIP_THRES))," fine-mapped SNPs (PIP>",PIP_THRES,") in ",length(snp_peaks <- unique(subset(trait_peaks_all,V5 > PIP_THRES)$V9))," peaks.\n\n"))
peak.df = subset(peak.data.df,peakid %in% snp_peaks)
fwrite(peak.df,paste0("/home/amarder/tmp/rbc.fm_peak.pip_",PIP_THRES,".txt"),quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)







