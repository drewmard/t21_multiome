# post_scent.R

library(data.table)
celltype_to_use="HSCs_all"
fDir="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/out_split"
fnums=paste0("/oak/stanford/groups/smontgom/amarder/t21_multiome/output/scent/input/info/",celltype_to_use,"_nums")
nums = fread(fnums,data.table = F,stringsAsFactors = F,header = F)

df.lst = list()
for (i in nums[,1]) {
  # print(i)
  f=paste0(fDir,"/",celltype_to_use,".",i,".txt")
  # if (!file.exists(f)) {print(i); next}
  df.lst[[i]] = fread(f,data.table = F,stringsAsFactors = F)
}
df = as.data.frame(do.call(rbind,df.lst))
df$fdr = p.adjust(df$pval,method='fdr')
sum(df$fdr < 0.2,na.rm=T)
fOut=paste0(fDir,"/all/",celltype_to_use,".all.txt")
fwrite(df,fOut,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)