# analyses use "TGV try again" data file, sheet named "redo"
# neonate size uses "S c1 neo size" sheet from neonate size excel file

library(survival)

library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(coxme)
library(dplyr)

redo$Parasite<-recode_factor(redo$Parasite, None="aNone")
redo$MomParasite<-recode_factor(redo$MomParasite, None="aNone")
redo$MomPXPara<-paste(redo$MomParasite,redo$Parasite)

####supfig
suplife<-coxph(Surv(deathday,edied)~ MomParasite*Parasite+strata(Clone), data=supfig)
summary(suplife)
test.ph<-cox.zph(life)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(life))


####EARLY MORTALITY
life<-coxph(Surv(deathday,edied)~ MomParasite*Parasite, data=redo)
summary(life)
test.ph<-cox.zph(life)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(life))

###confirming early death differences with a second strategy. restricted mean.
install.packages("survRM2")
library(survRM2)
#recode momparasite and parasite into new 0/1 variables with control as 0. micg as 1.
redo$armPara<-redo$Parasite
redo$armMomPara<-redo$MomParasite
redo$armPara<-recode_factor(redo$armPara, aNone="0")
redo$armPara<-recode_factor(redo$armPara, "O. pajunii"="1")
redo$armMomPara<-recode_factor(redo$armMomPara, aNone="0")
redo$armMomPara<-recode_factor(redo$armMomPara, MicG="1")

###restricted mean model. effect of Parasite
dd<-redo$deathday
d<-redo$died
a<-redo$armPara
rmst2(dd,d,a,tau=20, alpha=0.05)
rmst2(redo$deathday, redo$died, redo$armPara, tau=20, alpha=0.5)
rmst2(redo$deathday, redo$died, redo$armMomPara, tau=20, alpha=0.5)

###restricted mean model. effect of MomParasite
dd<-redo$deathday
d<-redo$died
a<-redo$armMomPara
rmst2(dd,d,a,tau=20, alpha=0.05)

Slifeplot<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomPXPara, data = redo),
                      xlab = "Days", 
                      ylab = "Survival Probability",
                      conf.int = FALSE, 
                      risk.table = FALSE,
                      censor=FALSE,
                      linetype = c("solid","solid","dashed","dashed"),
                      size = 1, 
                      xlim=c(0,60),
                      ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                      palette = c("BLACK", "GREY", "BLACK", "GREY"),
                      legend = c(.2,.1), 
                      legend.title = "Parasite Exposure",
                      legend.labs = c( "None", "None/MicG","MicG/None", "MicG/MicG"))


Slifeplot

###Address reviewer comment of whther this is just impact of mom immune activation by comparing to other pahtogens
###earliest death in metschmom group was day28, altered analysis to include deaths until day 28 to avoid infinite coeff.
####Early mortality across pathogens
immact$MomParasite<-recode_factor(immact$MomParasite, None="aNone")
#immact$MomParasite<-recode_factor(immact$MomParasite, amicg="MicG")
#immact$MomParasite<-recode_factor(immact$MomParasite, MicG="amicg")
immlife<-coxph(Surv(deathday,edied)~ MomParasite, data=immact)
summary(immlife)
test.ph<-cox.zph(immlife)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(life))

###confirming early death differences with a second strategy. restricted mean.
install.packages("survRM2")
library(survRM2)
#recode momparasite and parasite into new 0/1 variables with control as 0. micg as 1.
redo$armPara<-redo$Parasite
redo$armMomPara<-redo$MomParasite
redo$armPara<-recode_factor(redo$armPara, aNone="0")
redo$armPara<-recode_factor(redo$armPara, MicG="1")
redo$armMomPara<-recode_factor(redo$armMomPara, aNone="0")
redo$armMomPara<-recode_factor(redo$armMomPara, MicG="1")

###restricted mean model. effect of Parasite
dd<-redo$deathday
d<-redo$died
a<-redo$armPara
rmst2(dd,d,a,tau=20, alpha=0.05)

###restricted mean model. effect of MomParasite
dd<-redo$deathday
d<-redo$died
a<-redo$armMomPara
rmst2(dd,d,a,tau=20, alpha=0.05)

immlifeplot<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomParasite, data = immact),
                      xlab = "Days", 
                      ylab = "Survival Probability",
                      conf.int = FALSE, 
                      risk.table = FALSE,
                      censor=FALSE,
                      #linetype = c("solid","solid","dashed","dashed"),
                      size = 1, 
                      xlim=c(0,20),
                      ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                      palette = c("RED", "BLACK", "BLUE", "YELLOW"),
                      legend = c(.35,.15), 
                      legend.title = "Maternal Parasite Exposure",
                      legend.labs = c( "O. pajunii", "Control","P. ramosa", "M. bicuspidata"))


immlifeplot



#########################################################
#days from birth to first reproduction
timerepro<-coxph(Surv(C1day)~ MomParasite*Parasite, data=redo)
summary(timerepro)
test.ph<-cox.zph(timerepro)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(timerepro))



##########################################################
REPRODUCTION

###Create sum for total repro. inclusive of zeros
redo$C1forsum<-redo$C1
redo$C1forsum[is.na(redo$C1forsum)] <- 0
redo$C2forsum<-redo$C2
redo$C2forsum[is.na(redo$C2forsum)] <- 0
redo$C3forsum<-redo$C3
redo$C3forsum[is.na(redo$C3forsum)] <- 0
redo <- transform( redo, sumrepro = C1forsum + C2forsum +C3forsum)
redo$sumrepro[is.na(redo$sumrepro)] <- 0

#create sum for total repro. only animals which had 3 clutches
redo$C1for3c<-redo$C1
redo$C2for3c<-redo$C2
redo$C3for3c<-redo$C3
redo <- transform( redo, sum3c = C1for3c + C2for3c +C3for3c)
redo2 <- transform( redo,sum3c=sum3c )
redo2<-redo2[!is.na(redo2$sum3c),]

#Total repro (clutches 1-3). all animals
library(ARTool)
TR<-art(sumrepro ~MomParasite*Parasite, data = redo)
summary(TR)
anova(TR)
#art.con(TR, "MomParasite:Parasite") not needed, no interaction
aggregate(redo$sumrepro, list(redo$MomParasite), FUN=mean)
aggregate(redo$sumrepro, list(redo$MomParasite), FUN=se)

PanelB<-ggplot(data=redo, aes(x = MomParasite, y=sumrepro, fill = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values = colorRampPalette(c("WHITE","BLACK"))(2)) +
  labs(x = "Maternal Exposure", y = "S Early clutch neonate production") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="top")+ theme(text = element_text(size = 17))  
  

PanelB


#total repro. only animasl that had 3clutches
TR3<-art(sum3c ~MomParasite*Parasite, data = redo2)
summary(TR3)
anova(TR3)
aggregate(redo2$sum3c, list(redo2$MomParasite), FUN=mean)
aggregate(redo2$sum3c, list(redo2$MomParasite), FUN=se)
aggregate(redo2$sum3c, list(redo2$Parasite), FUN=mean)

###first clutch neonate abundance animals who had no first clutch = 0
c1w0<-art(C1forsum ~MomParasite*Parasite, data = redo)
summary(c1w0)
anova(c1w0)

###how many clutches did they have
manyclutches<-art(clutches ~MomParasite*Parasite, data = redo)
summary(manyclutches)
anova(manyclutches)
aggregate(redo$clutches, list(redo$MomParasite), FUN=mean)
aggregate(redo$clutches, list(redo$MomParasite), FUN=se)


PanelCa<-ggplot(data=redo, aes(x = MomParasite, y=clutches, fill = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values = colorRampPalette(c("WHITE","BLACK"))(2)) +
  labs(x = "Maternal Exposure", y = "S reproductive events") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="top")+ theme(text = element_text(size = 17))  


PanelCa

####size at first reproduction
redo1 <- transform( redo,c1size=c1size )
redo1<-redo1[!is.na(redo1$c1size),]
reprosize<-art(c1size ~MomParasite*Parasite, data = redo1)
summary(reprosize)
anova(reprosize)
aggregate(redo1$c1size, list(redo1$MomParasite), FUN=mean)
aggregate(redo1$c1size, list(redo1$MomParasite), FUN=sd)

#SE = SD/(sqrt(N))


#####neonate size
neoglm <- glm(Size ~ Gmapara*mpara + Gmapara*mpara/MomID,data=neosize)
summary(neoglm)
anova(neoglm, test="F")
##check residuals. 
##for variance= red lines should be centered on 0 or 1 and horiz.
##for normality qq should should align on diag.
par(mfrow=c(2,2))
plot(neoglm)
par(mfrow=c(1,1))
##normality. p<0.05=violation of assump.
shapiro.test(rstandard(neoglm))
###raw and transformed size does not fit well


library(ARTool)
neosize$MomID<-as.factor(neosize$MomID)
neosize$Gmapara<-as.factor(neosize$Gmapara)
neosize$mpara<-as.factor(neosize$mpara)
neo<-art(Size ~Gmapara*mpara, data = neosize)
summary(neo)
anova(neo)
art.con(neo, "Gmapara:mpara")
neosize$GmaXma<-paste(neosize$Gmapara,neosize$mpara)
aggregate(neosize$Size, list(neosize$GmaXma), FUN=mean)
aggregate(neosize$Size, list(neosize$GmaXma), FUN=se)

library("ggpubr")
ggboxplot(redo, x = "MomParasite", y = "clutches", color = "Parasite",
          palette = c("BLACK", "#1B9E77") )
library("ggpubr")
ggboxplot(redo1, x = "MomParasite", y = "c1size", color = "Parasite",
          palette = c("BLACK", "#1B9E77") )
library("ggpubr")
ggboxplot(neosize, x = "mpara", y = "Size", color = "Gmapara",
          palette = c("BLACK", "#1B9E77") )
library("ggpubr")
ggboxplot(redo, x = "MomParasite", y = "sumrepro", color = "Parasite",
          palette = c("BLACK", "#1B9E77") )




####time to infection
infectedday<-coxph(Surv(INFday,INFcensor)~ MomParasite, data=griddata)
summary(infectedday)
test.ph<-cox.zph(infectedday)
test.ph
par(mfrow=c(2,2))

##infection prevalence 
res <- prop.test(x = c(8, 4), n = c(8, 5))
res
##too few samples to get a reliable estimate

### infection intensity

sumdif <- wilcox.test(sum ~ MomParasite, data = griddata,
                      exact = FALSE)
sumdif
sumalive <- wilcox.test(scoreperdayalive ~ MomParasite, data = griddata, exact=FALSE)
sumalive
suminf <- wilcox.test(scoreperdayinfected ~ MomParasite, data = griddata,
                      exact = FALSE)
suminf


#####neonate weight
wt <- wilcox.test(indave ~ matpara, data = neowt,
                      exact = FALSE)
wt
library("ggpubr")
ggboxplot(neowt, x = "matpara", y = "indave", color = "matpara",
          palette = c("BLACK", "#1B9E77") )

t.test(neowt$indave ~ neowt$matpara, var.equal = TRUE)

###########################################
# B and D genotypes
###########################################
library(ARTool)

Bneosize$Gmapara<-as.factor(Bneosize$Gmapara)
Bneosize$mpara<-as.factor(Bneosize$mpara)
bneo<-art(Size ~Gmapara*mpara, data = Bneosize)
summary(bneo)
anova(bneo)
art.con(bneo, "Gmapara:mpara")
neosize$GmaXma<-paste(neosize$Gmapara,neosize$mpara)
aggregate(neosize$Size, list(neosize$GmaXma), FUN=mean)
aggregate(neosize$Size, list(neosize$GmaXma), FUN=se)
