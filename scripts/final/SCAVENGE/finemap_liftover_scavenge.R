library(data.table)
fileDir="/oak/stanford/groups/smontgom/amarder/t21_multiome/output/finemap/finemappedtraits_hg19"
dir.create(paste0(fileDir,"/hg19"),showWarnings = F)
dir.create(paste0(fileDir,"/hg38"),showWarnings = F)
fileName="ALL_FMdata.txt"
fileName="mono.PP001.bed"
flst = list.files(fileDir,pattern = ".txt|.bed",include.dirs = FALSE)
i = 2

options(scipen=999)
for (i in 1:length(flst)) {
  print(i)
  fileName = flst[i]
  fileNamePrefix=substring(fileName,1,nchar(fileName)-4)
  f=paste0(fileDir,"/",fileName)
  df = fread(f,data.table = F,stringsAsFactors = F)
  df = df[,c(1:5)]
  if (df[1,2]==df[1,3]) {
    df[,2] = df[,2]-1
  }
  f.out = paste0(fileDir,"/hg19/",fileNamePrefix,".hg19.bed")
  fwrite(df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = F)
}

for (i in 1:length(flst)) {
  print(i)
  fileName = flst[i]
  fileNamePrefix=substring(fileName,1,nchar(fileName)-4)
  f.hg19 = paste0(fileDir,"/hg19/",fileNamePrefix,".hg19.bed")
  f.hg38 = paste0(fileDir,"/hg38/",fileNamePrefix,".hg38.bed")
  f.unmapp = paste0(fileDir,"/hg38/",fileNamePrefix,".unmapp.bed")
  system(paste0('
echo "Running liftOver..."
HEADDIR=',fileDir,'
bedFile=',f.hg19,'
hg19ToHg38chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg19ToHg38.over.chain.gz
source /oak/stanford/groups/smontgom/amarder/bin/conda_init.sh
conda activate liftover
liftOver $bedFile $hg19ToHg38chain ',f.hg38,' ',f.unmapp,'
conda deactivate 
'))
}
options(scipen=0)

