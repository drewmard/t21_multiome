library(data.table)
library(ggplot2)
expr.data = list()
peak.data = list()
celltype_to_use = "HSCs_H"
fout = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out_info/expr.",celltype_to_use,".txt")
expr.data[[1]] = fread(fout,data.table = F,stringsAsFactors = F)
fout = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out_info/peak.",celltype_to_use,".txt")
peak.data[[1]] = fread(fout,data.table = F,stringsAsFactors = F)

celltype_to_use = "HSCs_T21"
fout = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out_info/expr.",celltype_to_use,".txt")
expr.data[[2]] = fread(fout,data.table = F,stringsAsFactors = F)
fout = paste0("/Users/andrewmarderstein/Documents/Research/t21_multiome/output/scent/out_info/peak.",celltype_to_use,".txt")
peak.data[[2]] = fread(fout,data.table = F,stringsAsFactors = F)

expr.data = as.data.frame(do.call(cbind,expr.data))
peak.data = as.data.frame(do.call(cbind,peak.data))

expr.data = expr.data[,c(1,2,4)]
peak.data = peak.data[,c(1,2,4)]
colnames(expr.data) = c("gene","expr_H","expr_t21")
colnames(peak.data) = c("peak","acc_H","acc_t21")

df = fread("~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs.DE.txt",data.table = F,stringsAsFactors = F)
df = merge(df,expr.data,by='gene')
df = merge(df,peak.data,by='peak')
df.sub = subset(df,fdr_H > 0.2 & fdr_t21_dn < 0.2)
# df.sub = subset(df,fdr_H < 0.2 & fdr_t21_dn > 0.2)
t.test(df.sub$acc_H,df.sub$acc_t21)
t.test(df.sub$expr_H,df.sub$expr_t21)
ggplot(df.sub,aes(x=acc_H,y=acc_t21)) + geom_point() + 
  geom_smooth(col='red') + theme_bw() + 
  geom_abline(slope = 1,intercept = 0,lty='dashed')
ggplot(df.sub,aes(x=expr_H,y=expr_t21)) + geom_point() + 
  geom_smooth(col='red') + theme_bw() + 
  geom_abline(slope = 1,intercept = 0,lty='dashed')

subset(df,fdr_H > 0.2 & fdr_t21 < 0.2 & P.Value_smATAC > 0.05)[1:2,]
ggplot(df,aes(x=-log10(boot_basic_p_t21),y=-log10(boot_basic_p_t21_dn))) + geom_point() + geom_smooth()
ggplot(df,aes(x=-log10(boot_basic_p_t21),y=-log10(boot_basic_p_H))) + geom_point() + geom_smooth()




