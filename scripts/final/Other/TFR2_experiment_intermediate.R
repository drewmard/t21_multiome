library(data.table)
library(ggplot2)
df1 = fread("~/Downloads/Ctrl_TFR2+.csv",data.table = F,stringsAsFactors = F)
df2 = fread("~/Downloads/TFR2_TFR2 OE .csv",data.table = F,stringsAsFactors = F)
df1$set="1"
df2$set="2"
tmp = as.data.frame(rbind(df1,df2))
ggplot(tmp,aes(fill=set,x=log10(`CD235`))) + geom_density(aes(fill=set),alpha=0.3)

median(df1$`CD235`);median(df2$`CD235`)
median(df1$`CD36`);median(df2$`CD36`)
median(df1$`CD71`);median(df2$`CD71`)
wilcox.test(df1$`CD235`,df2$`CD235`)$p.value
wilcox.test(df1$`CD36`,df2$`CD36`)$p.value
wilcox.test(df1$`CD71`,df2$`CD71`)$p.value
t.test(df1$`CD235`,df2$`CD235`)
t.test(df1$`CD36`,df2$`CD36`)
t.test(df1$`CD71`,df2$`CD71`)

mean(df2$`CD71`) / mean(df1$`CD71`)
mean(df2$`CD235`) / mean(df1$`CD235`)

t.test(df1$`CD235`,df2$`CD235`)
wilcox.test(df1$`CD235`,df2$`CD235`)$p.value
t.test(df1$`CD36`,df2$`CD36`)
wilcox.test(df1$`CD36`,df2$`CD36`)
wilcox.test(df1$`CD71`,df2$`CD71`)

# plot(log10(df1$`Comp-APC-A :: TFR2`),log10(df1$`CD235`))
# plot(log10(df1$`CD36`),log10(df1$`CD235`))
