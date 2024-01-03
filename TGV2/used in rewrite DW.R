#for supplement of TGV rewrite. DW clone analyses
# use TGV try again excel file, sheet redo D


library(survival)

library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(coxme)
library(dplyr)

redoD$Parasite<-recode_factor(redoD$Parasite, None="aNone")
redoD$MomParasite<-recode_factor(redoD$MomParasite, None="aNone")
redoD$MomPXPara<-paste(redoD$MomParasite,redoD$Parasite)

lifeD<-coxph(Surv(deathday,edied)~ MomParasite*Parasite, data=redoD)
summary(lifeD)
test.ph<-cox.zph(lifeD)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(lifeD))
#interaction infinite, cannot use this analysis. 

###checking early death differences with a second strategy. restricted mean.
install.packages("survRM2")
library(survRM2)
#recode momparasite and parasite into new 0/1 variables with control as 0. micg as 1.
redoB$armPara<-redo$Parasite
redoB$armMomPara<-redoB$MomParasite
redoB$armPara<-recode_factor(redoB$armPara, aNone="0")
redoB$armPara<-recode_factor(redoB$armPara, "O. pajunii"="1")
redoB$armMomPara<-recode_factor(redoB$armMomPara, aNone="0")
redoB$armMomPara<-recode_factor(redoB$armMomPara, MicG="1")
library(survRM2)
###restricted mean model. effect of Parasite
rmst2(redoD$deathday, redoD$died, redoD$armpara, tau=20, alpha=0.05)


###restricted mean model. effect of MomParasite
###redoD$arm = maternal pathogen treatment
rmst2(redoD$deathday, redoD$died, redoD$arm, tau=20, alpha=0.05)
###Indicates increased risk of early mortality due to maternal pathogen exposure

Dlifeplot<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomPXPara, data = redoD),
                      xlab = "Days", 
                      ylab = "Survival Probability",
                      conf.int = FALSE, 
                      risk.table = FALSE,
                      censor=FALSE,
                      linetype = c("solid","solid","dashed","dashed"),
                      size = 1, 
                      xlim=c(0,20),
                      ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                      palette = c("BLACK", "GREY", "BLACK", "GREY"),
                      legend = c(.2,.1), 
                      legend.title = "Parasite Exposure",
                      legend.labs = c( "None/None", "None/MicG","MicG/None", "MicG/MicG"))


Dlifeplot

Dlifeplotsimple<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomParasite, data = redoD),
                      xlab = "Days", 
                      ylab = "Survival Probability",
                      conf.int = FALSE, 
                      risk.table = FALSE,
                      censor=FALSE,
                      #linetype = c("solid","solid","dashed","dashed"),
                      size = 1, 
                      xlim=c(0,20),
                      ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                      #palette = c("BLACK", "GREY", "BLACK", "GREY"),
                      legend = c(.2,.1), 
                      legend.title = "Parasite Exposure")
                      #legend.labs = c( "None/None", "None/MicG","MicG/None", "MicG/MicG"))


Dlifeplotsimple

#########################################################
#days from birth to first reproduction
timereproD<-coxph(Surv(C1day)~ MomParasite*Parasite, data=redoD)
summary(timereproD)
test.ph<-cox.zph(timereproD)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(timereproD))




##########################################################
REPRODUCTION

###Create sum for total repro. inclusive of zeros
redoD$C1forsum<-redoD$C1
redoD$C1forsum[is.na(redoD$C1forsum)] <- 0
redoD$C2forsum<-redoD$C2
redoD$C2forsum[is.na(redoD$C2forsum)] <- 0
redoD$C3forsum<-redoD$C3
redoD$C3forsum[is.na(redoD$C3forsum)] <- 0
redoD <- transform( redoD, sumrepro = C1forsum + C2forsum +C3forsum)
redoD$sumrepro[is.na(redoD$sumrepro)] <- 0

#create sum for total repro. only animals which had 3 clutches
redoD$C1for3c<-redoD$C1
redoD$C2for3c<-redoD$C2
redoD$C3for3c<-redoD$C3
redoD <- transform( redoD, sum3c = C1for3c + C2for3c +C3for3c)
redo2D <- transform( redoD,sum3c=sum3c )
redo2D<-redo2D[!is.na(redo2D$sum3c),]


#USED IN SI. Total repro (clutches 1-3). all animals
library(ARTool)
TRD<-art(sumrepro ~MomParasite*Parasite, data = redoD)
summary(TRD)
anova(TRD)
#art.con(TR, "MomParasite:Parasite") not needed, no interaction
aggregate(redoD$sumrepro, list(redoD$Parasite), FUN=mean)
aggregate(redoD$sumrepro, list(redoD$Parasite), FUN=sd)

#total repro. only animasl that had 3clutches
TR3D<-art(sum3c ~MomParasite*Parasite, data = redo2D)
summary(TR3D)
anova(TR3D)
aggregate(redo2D$sum3c, list(redo2D$Parasite), FUN=mean)
aggregate(redo2D$sum3c, list(redo2D$MomParasite), FUN=se)
aggregate(redo2D$sum3c, list(redo2D$Parasite), FUN=mean)


#USED IN SI. repro only c1
redo21D<-redoD
redo21D$C1[is.na(redo21D$C1)] <- 0
TR3D<-art(C1 ~MomParasite*Parasite, data = redo21D)
summary(TR3D)
anova(TR3D)

####USED IN SI. size at first reproduction
redo1D <- transform( redoD,c1size=c1size )
redo1D<-redo1D[!is.na(redo1D$c1size),]
reprosizeD<-art(c1size ~MomParasite*Parasite, data = redo1D)
summary(reprosizeD)
anova(reprosizeD)
aggregate(redo1D$c1size, list(redo1D$MomParasite), FUN=mean)
aggregate(redo1D$c1size, list(redo1D$MomParasite), FUN=se)

###how many clutches did they have
redoD$MomParasite<-as.factor(redoD$MomParasite)
redoD$Parasite<-as.factor(redoD$Parasite)
manyclutchesD<-art(clutches ~MomParasite*Parasite, data = redoD)
summary(manyclutchesD)
anova(manyclutchesD)
art.con(manyclutchesD, "MomParasite:Parasite")
#under current exposure, maternal exposure decreases number of clutches

PanelCb<-ggplot(data=redoD, aes(x = MomParasite, y=clutches, fill = Parasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values = colorRampPalette(c("WHITE","BLACK"))(2)) +
  labs(x = "Maternal Exposure", y = "D reproductive events") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="top")+ theme(text = element_text(size = 17))  


PanelCb


aggregate(redoD$clutches, list(redoD$MomParasite), FUN=mean)
aggregate(redoD$clutches, list(redoD$MomParasite), FUN=se)
aggregate(redoD$clutches, list(redoD$MomPXPara), FUN=mean)
library(ggpubr)
ggboxplot(redoD, x = "MomParasite", y = "clutches", color = "Parasite",
          palette = c("#00AFBB", "#E7B800") )



