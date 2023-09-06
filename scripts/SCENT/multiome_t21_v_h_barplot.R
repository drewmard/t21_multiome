library(data.table)
x=data.frame(
  type=c("Healthy","Ts21","Ts21 (downsampled)"),
  n=c(191,1564,773)
)
x
library(ggplot2)
ggplot(x,aes(x=type,y=n,fill=type)) + geom_bar(stat='identity',col='black') + ggpubr::theme_pubr() +
  theme(axis.text.x = element_text(vjust=1,hjust=1,angle=30),
        axis.title.x = element_blank()) +
  scale_fill_manual(values=c("darkred","steelblue","steelblue4")) +
  guides(fill='none') +
  labs(y="# sig peak-gene links (FDR < 0.2)") +
  scale_y_continuous(breaks=seq(0,1600,by=200))

191-69

fisher.test(matrix(c(30000,191-69,773-69,69),2,2))$p.value
