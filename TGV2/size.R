

install.packages("ggplot2")
library(ggplot2)




####size at first reproduction
redo1 <- transform( redo,c1size=c1size )
redo1<-redo1[!is.na(redo1$c1size),]
reprosize<-art(c1size ~MomParasite*Parasite, data = redo1)
summary(reprosize)
anova(reprosize)
aggregate(redo1$c1size, list(redo1$MomParasite), FUN=mean)
aggregate(redo1$c1size, list(redo1$MomParasite), FUN=sd)
redo1
#SE = SD/(sqrt(N))

#####CODE FOR Size first repro figure

sz<-ggplot(redo1, aes(x=MomParasite, y=c1size  )) +
  geom_boxplot()+theme_classic()+ labs(title="",x="MomParasite", y = "c1size")
sz+scale_fill_manual(values=c("#636363", "#cccccc"))
