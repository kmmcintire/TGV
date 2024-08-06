##Analyses used in TGV manuscript. 8.5.2024
#title: Transgenerational pathogen effects: Maternal pathogen exposure reduces offspring fitness
#Submitted to Ecology for review


#EARLY MORTALITY - MAIN EXPERIMENT, early mortality - maternal Op vs current metsch
#TIME TO 1ST REPRO, NEONATE SIZE, # of clutches, size at 1st repro., burden
#REVIEWER SUGGESTED ADDITIONAL ANALYSES: early life maternal Op vs none/none. early life current Op vs current metsch 

# analyses use "tgv 8524" data file



library(survRM2)
library(dplyr)
library(survival)
library(lme4)
library(ggplot2)
library(ARTool)


###############################################################
#EARLY MORTALITY ANALYSES. reported in  supplemental stats table 8.2.24
#Main EXP sheet imported as "redoz"
#Analysis of 4 main groups: maternal parasite= Op or None. Current Parasite = Op or None
#recode momparasite and parasite into new 0/1 variables with control as 0. micg as 1.
redoz$armPara<-redoz$Parasite
redoz$armMomPara<-redoz$MomParasite
redoz$armPara<-recode_factor(redoz$armPara, None="0")
redoz$armPara<-recode_factor(redoz$armPara, "MicG"="1")
redoz$armMomPara<-recode_factor(redoz$armMomPara, None="0")
redoz$armMomPara<-recode_factor(redoz$armMomPara, MicG="1")

###restricted mean model. effect of Parasite. difference reported in table 8.2.24
rmst2(redoz$deathday, redoz$died, redoz$armPara, tau=20, alpha=0.5)

###restricted mean model. effect of Mom Parasite. difference reported in table 8.2.24
rmst2(redoz$deathday, redoz$died, redoz$armMomPara,,tau=20, alpha=0.05)


###################################
#Early mortality - TGV (maternal O.p./no current) VS no maternal/current METSCH
#reported in supplemental table. 8.2.24
#TGV vs Metsch sheet, imported as "tgvvmz"
# met=metsch infected. tgv=mom exposed to micg. 

#remove row 20= exposed but ultimately not infected with metsch
tormz <- tgvvmz[-c(20),]

tormz$arm<-tormz$tgvormet
tormz$arm<-recode_factor(tormz$arm, met="0")
tormz$arm<-recode_factor(tormz$arm, tgv="1")
###restricted mean model. effect of maternal micg vs current metsch. early life
rmst2(tormz$deathday, tormz$died, tormz$arm, tau=20, alpha=0.05)
###restricted mean model. effect of maternal micg vs current metsch. total life
rmst2(tormz$deathday, tormz$died, tormz$arm, tau=NULL, alpha=0.05)
#results indicate early life difference. no total life difference.

#########################################################
#days from birth to first reproduction. HR reported in supplement table 8.2.24
redoz$Parasite<-recode_factor(redoz$Parasite, None="aNone")
redoz$MomParasite<-recode_factor(redoz$MomParasite, None="aNone")
timereproz<-coxph(Surv(C1day)~ MomParasite*Parasite, data=redoz)
summary(timereproz)
test.ph<-cox.zph(timereproz)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(timereproz))
length(timereproz$residuals)


##########################################################
#REPRODUCTION - # of babies. in supplement table 8.2.24

###Create sum for total repro. inclusive of zeros
redoz$C1forsum<-redoz$C1
redoz$C1forsum[is.na(redoz$C1forsum)] <- 0
redoz$C2forsum<-redoz$C2
redoz$C2forsum[is.na(redoz$C2forsum)] <- 0
redoz$C3forsum<-redoz$C3
redoz$C3forsum[is.na(redoz$C3forsum)] <- 0
redoz <- transform( redoz, sumrepro = C1forsum + C2forsum +C3forsum)
redoz$sumrepro[is.na(redoz$sumrepro)] <- 0

#create sum for total repro. only animals which had 3 clutches
redoz$C1for3c<-redoz$C1
redoz$C2for3c<-redoz$C2
redoz$C3for3c<-redoz$C3
redoz <- transform( redoz, sum3c = C1for3c + C2for3c +C3for3c)
redo2z <- transform( redoz,sum3c=sum3c )
redo2z<-redo2z[!is.na(redo2z$sum3c),]

#Total repro (clutches 1-3). all animals
#in table 8.2.24
TRz<-art(sumrepro ~MomParasite*Parasite, data = redoz)
summary(TRz)
anova(TRz)
#art.con(TR, "MomParasite:Parasite") not needed, no interaction
aggregate(redoz$sumrepro, list(redoz$MomParasite), FUN=mean)
aggregate(redoz$sumrepro, list(redoz$MomParasite), FUN=se)
#total repro. only animasl that had 3clutches
#in table 8.2.24
TR3z<-art(sum3c ~MomParasite*Parasite, data = redo2z)
summary(TR3z)
anova(TR3z)
aggregate(redo2z$sum3c, list(redo2z$MomParasite), FUN=mean)
aggregate(redo2z$sum3c, list(redo2z$MomParasite), FUN=se)
aggregate(redo2z$sum3c, list(redo2z$Parasite), FUN=mean)
length(TR3z$residuals)


#############################################################
#NeoNate Size, artanova with mom as random grouping in supplement table, 8.2.24
neolmer <- lmer(Size ~ Gmapara*mpara+(1|MomID), data=NeosizeA)
summary(neolmer)
tdat <- data.frame(predicted=predict(neolmer), residual = residuals(neolmer), referrer=NeosizeA$MomID)
ggplot(tdat,aes(x=predicted,y=residual)) + geom_point() + geom_hline(yintercept=0, lty=3)
ggplot(tdat,aes(x=predicted,y=residual, colour=referrer)) + geom_point() + geom_hline(yintercept=0, lty=3) + theme(legend.position = "none")
ggplot(tdat,aes(x=residual)) + geom_histogram(bins=20, color="black")
ggplot(tdat,aes(sample=residual)) + stat_qq() + stat_qq_line()
#variance and normality issues, raw and transformed.
##opt for non parametric ART ANOVA
NeosizeA$MomID<-as.factor(NeosizeA$MomID)
NeosizeA$Gmapara<-as.factor(NeosizeA$Gmapara)
NeosizeA$mpara<-as.factor(NeosizeA$mpara)
NeosizeA$MomID<-as.factor(NeosizeA$MomID)
neo<-art(Size ~Gmapara*mpara +(1|MomID), data = NeosizeA)
summary(neo)
anova(neo)


###########################################
#infection - O.p burden - in supplemental table. 8.2.24
#griddata sheet imported as "griddataz"
#analysis of mean pathogen burden over the animal's life. 
#scoreperday = total observed parasite burden (sum of grid squares)/days of animal lived (birth to death)
sumalive <- wilcox.test(scoreperday ~ mpara, data = griddataz, exact=FALSE)
sumalive

################################################################
#reviewer Suggested analyses. in supplement table. 8.2.24
######
#EARLY MORTALITY- TGV vs CONTROL
#CREATE set with only no parasite animals
ctrlvtgvz<-redoz[redoz$Parasite=="aNone",]
ctrlvtgvz$arm<-ctrlvtgvz$MomParasite
ctrlvtgvz$arm<-recode_factor(ctrlvtgvz$arm, aNone="0")
ctrlvtgvz$arm<-recode_factor(ctrlvtgvz$arm, MicG="1")
###restricted mean model. effect of no exposure vs maternal micg. early life
rmst2(ctrlvtgvz$deathday, ctrlvtgvz$died, ctrlvtgvz$arm, tau=20, alpha=0.05)

#####
#EARLY MORTALITY - Current O.p vs Current Metsch
#create set
currentop<-redoz[redoz$MomParasite=="aNone",]
currentop<-currentop[currentop$Parasite=="MicG",]
currentop <- currentop[c("Parasite","deathday", "died")]
currentmetsch<-tormz[tormz$tgvormet=="met",] 
currentmetsch$Parasite<-currentmetsch$para
currentmetsch <- currentmetsch[c("Parasite","deathday", "died")]
currents<-rbind(currentop,currentmetsch)

currents$arm<-currents$Parasite
currents$arm<-recode_factor(currents$arm, "MicG"="0")
currents$arm<-recode_factor(currents$arm, Metsch="1")
###restricted mean model. effect of current micg vs current metsch. early life
rmst2(currents$deathday, currents$died, currents$arm, tau=20, alpha=0.05)
###restricted mean model. effect of current micg vs current metsch. total life
rmst2(currents$deathday, currents$died, currents$arm, tau=NULL, alpha=0.05)
