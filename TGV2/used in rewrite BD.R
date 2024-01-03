#for supplement of TGV rewrite. BD clone analyses
# use TGV try again excel file, sheet redo B


library(survival)

library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(coxme)
library(dplyr)

redoB$Parasite<-recode_factor(redoB$Parasite, None="aNone")
redoB$MomParasite<-recode_factor(redoB$MomParasite, None="aNone")
redoB$MomPXPara<-paste(redoB$MomParasite,redoB$Parasite)

lifeB<-coxph(Surv(deathday,died)~ MomParasite*Parasite, data=redoB)
summary(lifeB)
test.ph<-cox.zph(lifeB)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(lifeB))
####cannot include interaction as may be infinite. no death in 1 group.


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
time<-redoB$deathday
status<-redoB$died
armB<-redoB$armPara
status<-as.numeric(status)
time<-as.numeric(time)
armB<-as.numeric(armB)
rmst2(time, status,armB,tau=60, alpha=0.05)
rmst2(redoB$deathday, redoB$died, redoB$armPara, tau=20, alpha=0.05)


###restricted mean model. effect of MomParasite
###redoB$arm = maternal pathogen treatment
rmst2(redoB$deathday, redoB$died, redoB$arm, tau=20, alpha=0.05)

#########################################################
#days from birth to first reproduction
timereproB<-coxph(Surv(C1day)~ MomParasite*Parasite, data=redoB)
summary(timereproB)
test.ph<-cox.zph(timereproB)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(timereproB))

############
#TIME to infection
timeinfB<-coxph(Surv(INFday, INFcensor)~ mompara, data=gridB)
summary(timeinfB)
test.ph<-cox.zph(timeinfB)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(timereproB))


### infection intensity

sumdif <- wilcox.test(sum ~ mompara, data = gridB,
                      exact = FALSE)
sumdif
sumalive <- wilcox.test(scoreperdayalive ~ mompara, data = gridB, exact=FALSE)
sumalive
suminf <- wilcox.test(scoreperdayinfected ~ mompara, data = gridB,
                      exact = FALSE)
suminf


##########################################################
REPRODUCTION

###Create sum for total repro. inclusive of zeros
redoB$C1forsum<-redoB$C1
redoB$C1forsum[is.na(redoB$C1forsum)] <- 0
redoB$C2forsum<-redoB$C2
redoB$C2forsum[is.na(redoB$C2forsum)] <- 0
redoB$C3forsum<-redoB$C3
redoB$C3forsum[is.na(redoB$C3forsum)] <- 0
redoB <- transform( redoB, sumrepro = C1forsum + C2forsum +C3forsum)
redoB$sumrepro[is.na(redoB$sumrepro)] <- 0

#create sum for total repro. only animals which had 3 clutches
redoB$C1for3c<-redoB$C1
redoB$C2for3c<-redoB$C2
redoB$C3for3c<-redoB$C3
redoB <- transform( redoB, sum3c = C1for3c + C2for3c +C3for3c)
redo2B <- transform( redoB,sum3c=sum3c )
redo2B<-redo2B[!is.na(redo2B$sum3c),]


#USED IN SI. Total repro (clutches 1-3). all animals
library(ARTool)
TRB<-art(sumrepro ~MomParasite*Parasite, data = redoB)
summary(TRB)
anova(TRB)
#art.con(TR, "MomParasite:Parasite") not needed, no interaction
aggregate(redoB$sumrepro, list(redoB$Parasite), FUN=mean)
aggregate(redoB$sumrepro, list(redoB$Parasite), FUN=sd)

#total repro. only animasl that had 3clutches
TR3b<-art(sum3c ~MomParasite*Parasite, data = redo2B)
summary(TR3b)
anova(TR3b)
aggregate(redo2B$sum3c, list(redo2B$Parasite), FUN=mean)
aggregate(redo2B$sum3c, list(redo2B$MomParasite), FUN=se)
aggregate(redo2B$sum3c, list(redo2B$Parasite), FUN=mean)


#USED IN SI. repro only c1
redo21B<-redoB
redo21B$C1[is.na(redo21B$C1)] <- 0
TR3b<-art(C1 ~MomParasite*Parasite, data = redo21B)
summary(TR3b)
anova(TR3b)


####USED IN SI. size at first reproduction
redo1B <- transform( redoB,c1size=c1size )
redo1B<-redo1B[!is.na(redo1B$c1size),]
reprosize<-art(c1size ~MomParasite*Parasite, data = redo1B)
summary(reprosize)
anova(reprosize)
aggregate(redo1B$c1size, list(redo1B$MomParasite), FUN=mean)
aggregate(redo1B$c1size, list(redo1B$MomParasite), FUN=se)

###how many clutches did they have
redoB$MomParasite<-as.factor(redoB$MomParasite)
redoB$Parasite<-as.factor(redoB$Parasite)
manyclutchesB<-art(clutches ~MomParasite*Parasite, data = redoB)
summary(manyclutchesB)
anova(manyclutchesB)
aggregate(redo$clutches, list(redo$MomParasite), FUN=mean)
aggregate(redo$clutches, list(redo$MomParasite), FUN=se)