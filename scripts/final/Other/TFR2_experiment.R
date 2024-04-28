library(data.table)
library(ggplot2)
df1 = fread("~/Downloads/Ctrl_TFR2 +.csv",data.table = F,stringsAsFactors = F)
df2 = fread("~/Downloads/TFR2_OE__TFR2 ++ .csv",data.table = F,stringsAsFactors = F)
df1$set="1"
df2$set="2"
tmp = as.data.frame(rbind(df1,df2))
ggplot(tmp,aes(fill=set,x=log10(`Comp-Brilliant UV 395-A :: CD235`))) + geom_density(aes(fill=set),alpha=0.3)

median(df1$`Comp-Brilliant UV 395-A :: CD235`);median(df2$`Comp-Brilliant UV 395-A :: CD235`)
median(df1$`Comp-PE-A :: CD36`);median(df2$`Comp-PE-A :: CD36`)
median(df1$`Comp-BV785-A :: CD71`);median(df2$`Comp-BV785-A :: CD71`)
wilcox.test(df1$`Comp-Brilliant UV 395-A :: CD235`,df2$`Comp-Brilliant UV 395-A :: CD235`)$p.value
wilcox.test(df1$`Comp-PE-A :: CD36`,df2$`Comp-PE-A :: CD36`)$p.value
wilcox.test(df1$`Comp-BV785-A :: CD71`,df2$`Comp-BV785-A :: CD71`)$p.value

t.test(df1$`Comp-Brilliant UV 395-A :: CD235`,df2$`Comp-Brilliant UV 395-A :: CD235`)
wilcox.test(df1$`Comp-Brilliant UV 395-A :: CD235`,df2$`Comp-Brilliant UV 395-A :: CD235`)$p.value
t.test(df1$`Comp-PE-A :: CD36`,df2$`Comp-PE-A :: CD36`)
wilcox.test(df1$`Comp-PE-A :: CD36`,df2$`Comp-PE-A :: CD36`)
wilcox.test(df1$`Comp-BV785-A :: CD71`,df2$`Comp-BV785-A :: CD71`)

plot(log10(df1$`Comp-APC-A :: TFR2`),log10(df1$`Comp-Brilliant UV 395-A :: CD235`))
plot(log10(df1$`Comp-PE-A :: CD36`),log10(df1$`Comp-Brilliant UV 395-A :: CD235`))
