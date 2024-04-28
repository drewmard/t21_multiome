library(data.table)

# 0. SCENT!
celltype_to_use="HSCs_H"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out_split/all/",celltype_to_use,".all.txt")
res.df.h = fread(output_file,data.table = F,stringsAsFactors = F)
colnames(res.df.h)[4:9] = paste0(colnames(res.df.h)[4:9],"_H")
res.df.h = res.df.h[,2:9]

celltype_to_use="HSCs_T21"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out_split/all/",celltype_to_use,".all.txt")
res.df.t21 = fread(output_file,data.table = F,stringsAsFactors = F)
colnames(res.df.t21)[4:9] = paste0(colnames(res.df.t21)[4:9],"_t21")
res.df.t21 = res.df.t21[,2:9]

celltype_to_use="HSCs_T21.down"
output_file = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out_split/all/",celltype_to_use,".all.txt")
res.df.t21.down = fread(output_file,data.table = F,stringsAsFactors = F)
colnames(res.df.t21.down)[4:9] = paste0(colnames(res.df.t21.down)[4:9],"_t21_dn")
res.df.t21.down = res.df.t21.down[,2:9]

res.df.mg = merge(res.df.h,res.df.t21,all=TRUE,by=c("gene","peak"))
res.df.mg = merge(res.df.mg,res.df.t21.down,all=TRUE,by=c("gene","peak"))

# 1. find trait-associated SNPs 
traitName="rbc"
trait_fm = fread(paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out1/",traitName,"_peaks_overlap.bed"),data.table = F,stringsAsFactors = F)
trait_fm.sub = trait_fm[,c("V1","V3","V4","V5","V9")]
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")
trait_fm.sub = subset(trait_fm,V5 > 0.2)[,c("V1","V3","V4","V5","V9")]
dim(trait_fm.sub)
nrow(subset(trait_fm,V5 > 0.9)[,c("V1","V3","V4","V5","V9")])
colnames(trait_fm.sub) = c("chr","snp_pos","rsid","pip","peak")

# 2. find the subset of peaks that are linked to genes
df.sub = merge(trait_fm.sub,res.df.mg,by="peak",all.y=TRUE)

# 3. 
small_ATAC_pb = fread("~/Downloads/pb_de_atac.hsc.txt",data.table = F,stringsAsFactors = F)
large_RNA_pb = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/Liver.HSCs_MPPs.sample.txt",data.table = F,stringsAsFactors = F)
small_RNA_pb = fread("~/Downloads/pb_de.hsc.txt",data.table = F,stringsAsFactors = F)
cyc_HSC_pb = fread("~/Documents/Research/t21-proj/out/full/DE_pb_leiden_names/HSC_v_CyclingHSC.txt",data.table = F,stringsAsFactors = F)
cyc_HSC_pb$P.Value.rank = rank(cyc_HSC_pb$P.Value)/nrow(cyc_HSC_pb)
colnames(small_ATAC_pb) = paste0(colnames(small_ATAC_pb),"_smATAC")
colnames(large_RNA_pb) = paste0(colnames(large_RNA_pb),"_lgRNA")
colnames(small_RNA_pb) = paste0(colnames(small_RNA_pb),"_smRNA")
colnames(cyc_HSC_pb) = paste0(colnames(cyc_HSC_pb),"_cycHSC")
df.sub2 = merge(df.sub,small_ATAC_pb,by.x='peak',by.y='names_smATAC',all.x=T)
df.sub2 = merge(df.sub2,small_RNA_pb,by.x='gene',by.y='names_smRNA',all.x=T)
df.sub2 = merge(df.sub2,large_RNA_pb,by.x='gene',by.y='names_lgRNA',all.x=T)
df.sub2 = merge(df.sub2,cyc_HSC_pb,by.x='gene',by.y='names_cycHSC',all.x=T)
df.sub2$fm = !(df.sub2$pip <= 0.2 | is.na(df.sub2$pip))

# to remove duplicates...
df.sub2 = df.sub2[order(-df.sub2$pip,df.sub2$p_t21),]
df.sub2 = df.sub2[!duplicated(df.sub2[,c("peak","gene")]),]
df.sub2 = df.sub2[order(df.sub2$gene,df.sub2$peak),]

fwrite(df.sub2,"~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs.DE.txt",quote = F,na = "NA",sep = "\t",row.names = F,col.names = T)

t.test(-log10(df.sub2$boot_basic_p_t21),-log10(df.sub2$boot_basic_p_H))
t.test(abs(df.sub2$beta_H),abs(df.sub2$beta_t21))
tmp = subset(df.sub2,!is.na(beta_H) & !is.na(beta_t21))
t.test(abs(tmp$beta_H),abs(tmp$beta_t21))


