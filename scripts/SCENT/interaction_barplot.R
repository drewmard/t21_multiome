library(ggplot2)
tmp = data.frame(
  Disease=c("Disomy","Disomy","Ts21","Ts21"),
  Interaction=c("Sig.","Not Sig.","Sig.","Not Sig."),
  Proportion=c(0.621,1-0.621,0.024,1-0.024)
)
g=ggplot(tmp,aes(x=Disease,fill=Interaction,y=Proportion)) + geom_bar(stat='identity',position = 'dodge',col='black') + ggpubr::theme_pubr() +
  scale_fill_manual(values=c("darkgrey","darkred"))
f.out="~/Downloads/interaction_plot.pdf"
pdf(f.out,width=4,height=3.3)
print(g)
dev.off()