library(data.table)

# celltype_to_use = "HSCs_T21"
# projName="tmparm"

assemble_parallel_files <- function(projName,celltype_to_use,read_files=FALSE) {
  
  
  if (read_files) {
    fDir = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input"
    fDirout = "/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/chunkinfo_split"
    f.chunkinfo = paste0(fDir,"/chunkinfo/",celltype_to_use,".chunkinfo.txt")
    f.atac_out = paste0(fDir,"/atac/",celltype_to_use,".atac.rds")
    f.rna_out = paste0(fDir,"/rna/",celltype_to_use,".rna.rds")
    
    chunkinfo = fread(f.chunkinfo,data.table = F,stringsAsFactors = F)
    atac.all = readRDS(f.atac_out)
    mrna = readRDS(f.rna_out)
  }
  
  tmpDir = paste0("/oak/stanford/groups/smontgom/amarder/tmp/",projName)
  unlink(paste0(tmpDir),recursive = TRUE)
  dir.create(paste0(tmpDir,"/"),showWarnings = FALSE)
  dir.create(paste0(tmpDir,"/atac"),showWarnings = FALSE)
  dir.create(paste0(tmpDir,"/rna"),showWarnings = FALSE)
  dir.create(paste0(tmpDir,"/chunkinfo_split"),showWarnings = FALSE)
  dir.create(paste0(fDir,"/info"),showWarnings = FALSE)
  if (file.exists(paste0(fDir,"/info/nums"))) {file.remove(paste0(fDir,"/info/nums"),showWarnings = FALSE)}
  
  rng = seq(1,nrow(chunkinfo),by=500)
  fwrite(data.frame(1:length(rng)),paste0(fDir,"/info/",celltype_to_use,"_nums"),row.names = F,col.names = F,quote = F,na = "NA",sep = "\t")
  
  for (i in 1:length(rng)) {
    start = rng[i]
    if (i==length(rng)) {
      end = nrow(chunkinfo)
    } else {
      end = rng[i+1] - 1
    }
    
    print(paste(start,end))
    
    chunkinfo.split = chunkinfo[start:end,]
    # fout.chunkinfo = paste0(fDir,"/chunkinfo_split/",celltype_to_use,".chunkinfo.",i,".txt")
    # fwrite(chunkinfo.split,fout.chunkinfo,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
    
    atac.sub = atac.all[unique(chunkinfo.split$peak),]
    mrna.sub = mrna[unique(chunkinfo.split$gene),]
    
    fout.atac = paste0(tmpDir,"/atac/",celltype_to_use,".atac.",i,".rds")
    fout.rna = paste0(tmpDir,"/rna/",celltype_to_use,".rna.",i,".rds")
    fout.chunkinfo = paste0(tmpDir,"/chunkinfo_split/",celltype_to_use,".chunkinfo.",i,".txt")
    
    saveRDS(atac.sub,fout.atac)
    saveRDS(mrna.sub,fout.rna)
    fwrite(chunkinfo.split,fout.chunkinfo,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  }
  
}

# Rscript scent_v3.R $celltype_to_use $i
# sbatch scent_v3.R

# --account=smontgom --partition=batch --time=4-1:00:00 --mem=32G --nodes=16

