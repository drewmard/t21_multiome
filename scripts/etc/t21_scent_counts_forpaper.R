library("exact2x2")
library(data.table)

df = fread("~/Documents/Research/t21_multiome/output/scent/out_split/all/HSCs.DE.int.txt",data.table = F,stringsAsFactors = F)
nrow(subset(df,fdr_t21 < 0.2))
nrow(subset(df,fdr_t21_dn < 0.2))
nrow(subset(df,fdr_H < 0.2))
nrow(subset(df,fdr_t21 < 0.2 & fdr_H < 0.2))
table(subset(df,fdr_t21 < 0.2 & fdr_H > 0.2)$fdr_int < 0.2)
table((subset(df,fdr_H < 0.2 & fdr_t21 > 0.2))$fdr_int < 0.2)
table((subset(df,fdr_H < 0.2 & fdr_t21 < 0.2))$fdr_int < 0.2)
table((subset(df,fdr_H < 0.2 & fdr_t21 < 0.2 & sign(beta_H)!=sign(beta_t21)))$fdr_int < 0.2)
(subset(df,fdr_H < 0.2 & fdr_t21 < 0.2 & sign(beta_H)!=sign(beta_t21)))

tmp = subset(df,fdr_t21 < 0.2 & fdr_H > 0.2 & fdr_int > 0.2)
mean(tmp$adj.P.Val_smATAC < 0.1)
mean(tmp$adj.P.Val_lgRNA < 0.1,na.rm = T)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H > 0.2)
mean(tmp$adj.P.Val_smATAC < 0.1)
mean(tmp$adj.P.Val_lgRNA < 0.1,na.rm = T)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H < 0.2 & fdr_int > 0.2)
mean(tmp$adj.P.Val_smATAC < 0.1)
mean(tmp$adj.P.Val_lgRNA < 0.1,na.rm = T)
table(tmp$adj.P.Val_lgRNA < 0.1)

tmp = subset(df,fdr_t21 < 0.2 & fdr_H > 0.2 & fdr_int > 0.2)
f0=mean(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0,na.rm = T);f0
y=table(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H > 0.2)
f=mean(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0,na.rm = T);f
f0/f;binom.test(y[2],sum(y),f)

tmp = subset(df,fdr_t21 < 0.2 & fdr_H > 0.2 & fdr_int > 0.2)
f0=mean(tmp$P.Value_smATAC < 0.05 & tmp$logFC_smATAC > 0,na.rm = T);f0
y=table(tmp$P.Value_smATAC < 0.05 & tmp$logFC_smATAC > 0)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H > 0.2)
f=mean(tmp$P.Value_smATAC < 0.05 & tmp$logFC_smATAC > 0,na.rm = T);f
f0/f;binom.test(y[2],sum(y),f)

tmp = subset(df,fdr_t21 > 0.2 & fdr_H < 0.2 & fdr_int > 0.2)
f0=mean(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0,na.rm = T);f0
y=table(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H > 0.2)
f=mean(tmp$P.Value_lgRNA < 0.05 & tmp$logFC_lgRNA > 0,na.rm = T);f
f0/f;binom.test(y[2],sum(y),f)

tmp = subset(df,fdr_t21 > 0.2 & fdr_H < 0.2 & fdr_int > 0.2)
f0=mean(tmp$P.Value_smATAC < 0.05 & tmp$logFC_smATAC > 0,na.rm = T);f0
y=table(tmp$P.Value_smATAC < 0.05 & tmp$logFC_smATAC > 0)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H > 0.2)
f=mean(tmp$P.Value_smATAC < 0.05 & tmp$logFC_smATAC > 0,na.rm = T);f
f0/f;binom.test(y[2],sum(y),f)




mean(tmp$adj.P.Val_lgRNA < 0.1,na.rm = T)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H > 0.2)
mean(tmp$adj.P.Val_smATAC < 0.1)
mean(tmp$adj.P.Val_lgRNA < 0.1,na.rm = T)
tmp = subset(df,fdr_t21 > 0.2 & fdr_H < 0.2 & fdr_int > 0.2)
mean(tmp$adj.P.Val_smATAC < 0.1)
mean(tmp$adj.P.Val_lgRNA < 0.1,na.rm = T)
table(tmp$adj.P.Val_lgRNA < 0.1)


tmp = subset(df,fdr_t21 > 0.2 & fdr_H > 0.2)
mean(tmp$adj.P.Val_smATAC < 0.2)

table(p.adjust((subset(df,fdr_H < 0.2 & fdr_t21 > 0.2))$p_int,method = 'fdr') < 0.2)
table(p.adjust((subset(df,fdr_H > 0.2 & fdr_t21 < 0.2))$p_int,method = 'fdr') < 0.2)


nrow(subset(df,fdr_H < 0.2 & fdr_t21 > 0.2))
nrow(subset(df,fdr_H < 0.2 & fdr_t21 > 0.2))

table((subset(df,fdr_H > 0.2 & fdr_t21 < 0.2))$fdr_int < 0.2)


table((subset(df,fdr_H < 0.2 & fdr_t21 > 0.2))$p_int < 0.05)
table((subset(df,fdr_H < 0.2 & fdr_t21 > 0.2))$p_int < 0.05)
table((subset(df,fdr_H > 0.2 & fdr_t21 < 0.2))$p_int < 0.05)
table((subset(df,fdr_H > 0.2 & fdr_t21 < 0.2))$boot_basic_p_int < 0.05)
table((subset(df,fdr_H > 0.2 & fdr_t21 < 0.2))$fdr_int < 0.2)

table(subset(df,fdr_h < 0.2 & fdr_t21 > 0.2)$fdr_int < 0.2)

df.sub = subset(df,fdr_H < 0.2 & fdr_int < 0.2)

nrow(df.sub)

