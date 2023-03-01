###TGV MANU

####get controlbabies dataset. 2/28
install.packages("MASS")
library(MASS)
#get the dataset
library(readxl)
controlbabies <- read_excel("Data/transgen data.xlsx", 
                            sheet = "control babies", col_types = c("numeric", 
                                                                    "text", "text", "text", "skip", "text", 
                                                                    "text", "skip", "numeric", "numeric", 
                                                                    "skip", "numeric", "numeric", "skip", 
                                                                    "numeric", "numeric", "skip", "numeric", 
                                                                    "skip", "numeric", "numeric", "skip"))
View(controlbabies)

###TOTAL REPRODUCTION. 2/28
#for first 3 clutches
#inclusive of all animals, even those with 0. 
#make new sumrepro variable, which sums the first 3 clutches. 
#by replacing NA with 0 for clutch variables and summing
controlbabies$C1forsum<-controlbabies$C1
controlbabies$C1forsum[is.na(controlbabies$C1forsum)] <- 0
controlbabies$C2forsum<-controlbabies$C2
controlbabies$C2forsum[is.na(controlbabies$C2forsum)] <- 0
controlbabies$C3forsum<-controlbabies$C3
controlbabies$C3forsum[is.na(controlbabies$C3forsum)] <- 0
controlbabies <- transform( controlbabies, sumrepro = C1forsum + C2forsum +C3forsum)
controlbabies$sumrepro[is.na(controlbabies$sumrepro)] <- 0
controlbabies <- transform( controlbabies, C1C2repro = C1forsum + C2forsum)
controlbabies$C1C2repro[is.na(controlbabies$C1C2repro)] <- 0
###FECUNDITY, GLOBAL 
#glm poisson
totalrepropois <- glm(sumrepro ~ MomParasite *Clone, data = controlbabies,  family = "poisson")
summary(totalrepropois)
##reside deviance much larger than DF, variance/mean>1, all signs point to overdispersion
###should not use poisson. should try NB
##GLM for total reproduction of all animals. NB used on count data.
##find theta for the nb analysis
fit = glm.nb(sumrepro ~ 1, data=controlbabies); fit$theta
#this above will return the theta value to be put in the models below
totrepronb<-glm.nb(sumrepro~ MomParasite * Clone , data=controlbabies, init.theta=1.757038, link=log)
summary(totrepronb)
par(mfrow = c(2,2))
plot(totrepronb)
anova(totrepronb, test = "Rao")
library("lmtest")
lrtest(totrepronb, totalrepropois)
###logliklihood of Negbin greater than poisson
###CONCLUSION=== nb reported in tgv manu supplemental table. no sig effect of maternal exposure or clone on total reproduction



####FECUNDITY FOR ONLY STD, 2/28
#includes animals producing 0
#create sum of first three clutches
transgen_data$C1forsum<-transgen_data$C1
transgen_data$C1forsum[is.na(transgen_data$C1forsum)] <- 0
transgen_data$C2forsum<-transgen_data$C2
transgen_data$C2forsum[is.na(transgen_data$C2forsum)] <- 0
transgen_data$C3forsum<-transgen_data$C3
transgen_data$C3forsum[is.na(transgen_data$C3forsum)] <- 0
transgen_data<- transform( transgen_data, sumrepro = C1forsum + C2forsum +C3forsum)
transgen_data$sumrepro[is.na(transgen_data$sumrepro)] <- 0
#mean total by maternal parasite
aggregate(transgen_data$sumrepro, list(transgen_data$MomParasite), FUN=mean) 
#neg bin glm for STD clone total repro
library(MASS)
totreprostd<-glm.nb(sumrepro~ MomParasite, data=transgen_data)
summary(totreprostd)
par(mfrow = c(2,2))
plot (totreprostd)
qqnorm(resid(totreprostd))
qqline(resid(totreprostd))
anova(totreprostd, test="Rao")
##suggests poor fit
##OPT FOR NONPARAMETRIC. KRUSKAL. 2/28
##REPORTED IN SUPP. TABLE
kruskal.test(sumrepro~ MomParasite, data=transgen_data)


####C1 NEONATE ABUNDANCE, 2/28
#glm for neonate abundance produced in first clutch
#includes animals producing 0
#MomParasite=parasite treatment of the mother of the daphnia creating this clutch
C1repro<-glm(C1forsum~ MomParasite * Clone , data=controlbabies)
anova(C1repro, test="F")
par(mfrow = c(2,2))
plot (C1repro)
qqnorm(resid(C1repro))
qqline(resid(C1repro))
#GLM normality questionable. try nonparametric analysis
# run the ART procedure for nonparametric analysis of C1neonate abundance
##ART REPORTED IN SUPPLEMENTAL TABLE 2/28
library(ARTool)
C1tot <- read.csv("C:\\Users\\krism\\Desktop\\transgenerational virulence\\transgen virulence\\Data\\controlbabiecsv.csv") # assumes file is in working directory
C1tot$C1forsum<-C1tot$C1
C1tot$C1forsum[is.na(C1tot$C1forsum)] <- 0
C1number = art(C1 ~ factor(Mpara) * factor(Clone) , data=C1tot) # linear mixed model syntax; see lme4::lmer
anova(C1number)
##non parametric reported in TGV manu


#####SIZE AT FIRST REPRODUCTION, 2/28
#GLMM  for mother's size at C1. Size log transformed
##REPORTED IN SUPPLEMENTAL TABLE
c1mm <- glm(log(c1size) ~ MomParasite * Clone, data = controlbabies)
summary(c1mm)
anova(c1size, test="F")
par(mfrow = c(2,2))
plot (c1mm)
qqnorm(resid(c1mm))
qqline(resid(c1mm))


###NEONATE SIZE, 2/28
##fit size distribution
library(fitdistrplus)
descdist(data=babysize$Size, discrete=FALSE)
descdist(data=babysize$Size, discrete=FALSE, boot=1000)
#descdist overlaps most with beta, but cannot be beta, also decent match for normal dist.
sizeranX <- lmer(log(Size)~GMApara * Clone + (1|MomID:Clone), data=babysize)
#test for main effect of GMApara
anova(sizeranX, test="Chisq")
neosize<-glm(log(Size)~GMApara*Clone, data=babysize)
anova(neosize, test="F")

##MICG TOTAL SPORE PRODUCTION, BURDEN. 2/28 
#analysis of total grid score
#chose nonparametric analysis, because resids did not conform well to simple distr.
# download and install the ARTool package
install.packages("ARTool")
# load the ARTool library into memory
library(ARTool)
# read a data table into variable 'df'
spores <- read.csv("C:\\Users\\krism\\Desktop\\transgenerational virulence\\transgen virulence\\Data\\micgspores.csv") # assumes file is in working directory
gridsum = art(sum ~ factor(mpara) * factor(Clone) , data=spores) # linear mixed model syntax; see lme4::lmer
summary(gridsum)
#Ftest not precisely zero, but approximately zero. acceptable.
anova(gridsum)
##gridsum reported in TGV
gridspd = art(scoreperday ~ factor(mpara) * factor(Clone) , data=spores) # linear mixed model syntax; see lme4::lmer
summary(gridspd)
###Ftests not close to zero, do not use
gridspid = art(scoreperinfectedday ~ factor(mpara) * factor(Clone) , data=spores) # linear mixed model syntax; see lme4::lmer
summary(gridspid)
# Ftests not close to zero, do not use


####INFECTION PREVALENCE, 2/28
####comparing models with and without interaction. not better with interaction.
###inf1tot used in tgv manuscript
micgbabies$INFtot<-micgbabies$INF
micgbabies$INFtot[is.na(micgbabies$INFtot)] <- 0
inf1tot <- glm(INFtot~ MomParasite * Clone, family = binomial, data = micgbabies)
summary(inf1tot)
par(mfrow = c(2,2))
plot(inf1tot)
inf2tot <- glm(INFtot ~ MomParasite + Clone, family = binomial, data = micgbabies)
summary(inf2tot)
anova(inf1tot, inf2tot, test = "Chisq")
anova(inf1tot, test = "Rao")
#ftest cannot be used with binomial, rao an alternative to wald



#####SURVIVAL- EARLY MORTALITY

#install and load packages for cox prop hazard and K-M
install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("cli")
install.packages("tidyverse")
install.packages("gridExtra")
library(tidyverse)
library(cli)
library(ggplot2)
library(survival)
library(ggpubr)
library(dplyr)
library(survminer)
library(gridExtra)




#get the dataset FOR MORTALITY ANALYSIS
library(readxl)
controlbabies <- read_excel("Data/transgen data.xlsx", 
                            sheet = "control babies", col_types = c("numeric", 
                                                                    "text", "text", "text", "skip", "text", 
                                                                    "text", "skip", "numeric", "numeric", 
                                                                    "skip", "numeric", "numeric", "skip", 
                                                                    "numeric", "numeric", "skip", "numeric", 
                                                                    "skip", "numeric", "numeric", "skip"))
View(controlbabies)
###TIME TO DEATH FOR EARLY LIFE 2/28
# cox proportional hazard for Control daughters, comparing control and micg moms, with clone included seperately
#IN S1TABLE
cphbymompclone<-coxph(Surv(deathday, censor)~Mpara + Clone, data=controlbabies, method="breslow")
summary(cphbymompclone)
#test for assumptions of prop hazards. p>0.05 means no departure from assumptions. 
test.ph <- cox.zph(cphbymompclone)
test.ph
#results indicate model fulfills assumption of prop. hazards


#####TIME to FIRST REPRODUCTION & INFECTION. 2/28
#COX prop hazard for time to reproduction
#IN S1 TABLE
cphforrepro<-coxph(Surv(C1day)~Mpara + strata(Clone) , data=controlbabies)
summary(cphforrepro)
test.ph<-cox.zph(cphforrepro)
test.ph
##COX prop hazard for time to infection
#IN S1 TABLE
timetoinf<-coxph(Surv(INFday,INFcensor)~Mpara , data=micgbabies)
summary(timetoinf)
test.ph <- cox.zph(timetoinf)
test.ph



##tgvlife=dataset with control mom and micgmom babies
##explife=dataset with controlmom and micGmom babies that have been exposed
explife<-tgvlife[tgvlife$Parasite=="MicG",]
stdlife<-explife[explife$Clone=="STD",]
dwlife<-explife[explife$Clone=="DW",]
bdlife<-explife[explife$Clone=="BD",]
complife<-tgvlife[tgvlife$MomParasite=="MicG",]


expcox<-coxph(Surv(deathday,censor)~MomParasite*strata(Clone), data=explife)
summary(expcox)
test.ph <- cox.zph(expcox)
test.ph
expplot<-ggsurvplot(fit=survfit(Surv(deathday,Ecensor)~MomParasite, data =stdlife),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = TRUE, 
                   risk.table = FALSE,
                   #linetype = c("solid","solid","solid", "dotted", "dotted", "dotted"),
                   size = 1, 
                   xlim=c(0,20),
                   #palette = c("#5e3c99", 
                           #   "#1B9E77", "#e66101", "#5e3c99", 
                            #  "#1B9E77", "#e66101"), 
                   legend = c(.15,.2),
                   legend.title = "STDE MomP")
                   #legend.labs = c("MicG BD", 
                               #  "MicG DW","MicG STD", "None BD", "None DW", "None STD"))
expplot

survdiff(Surv(deathday, censor) ~ MomParasite, data = explife)
survfit2(Surv(deathday, censor) ~ MomParasite, data = stdlife) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )
library(ggsurvfit)
install.packages("ggsurvfit")
install.packages("Greg")
library(Greg)
library(dplyr)
regular_model <- coxph(Surv(deathday, censor) ~
                         MomParasite*Clone,
                       data = explife)
summary(regular_model)

  timeSplitter(data= explife, by = 20,
               event_var = "censor",
               
               time_var = "deathday")
              # time_related_vars = c("MomParasite", "Clone"))

interval_model <-
  update(regular_model,
         Surv(Start_time, Stop_time, censor) ~ .,
         data = explife)

summary(interval_model)



