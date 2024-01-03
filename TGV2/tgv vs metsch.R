#for supplement of TGV rewrite, effect of TGV vs effect of known virulent early killing pathogen
#maternal micg vs current metsch


library(survival)

library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(coxme)
library(dplyr)
#tgvormet variable. met=metsch infected. tgv=mom exposed to micg. anone=no maternal or current exposure

#remove row 20 as this animal was exposed but ultimately not infected with metsch
# remove rows in r by row number
tgvvm2 <- tgvvm[-c(20),]
tgvvm2 <- tgvvm2[-c(18),]

#create new combined variable
tgvvm$MomPXPara<-paste(tgvvm$mompath,tgvvm$para)


lifesbad<-coxph(Surv(deathday,died)~ tgvormet, data=tgvvm2)
summary(lifesbad)
test.ph<-cox.zph(lifesbad)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(lifesbad))
#test.ph indicates not proportional hazard, therefore analysis results are unreliable
#opt for restricted mean time survival analysis below


###checking early death differences with a second strategy. restricted mean.
#this analysis will require 3 pairwise comparisons. tgv vs met. met vs none. tgv vs none.
install.packages("survRM2")
library(survRM2)

#recode tgvormet into new 0/1 variables with metsch as 0. maternalmicg as 1.
tgvvm2$arm<-tgvvm2$tgvormet
tgvvm2$arm<-recode_factor(tgvvm2$arm, met="0")
tgvvm2$arm<-recode_factor(tgvvm2$arm, anone="1")
###restricted mean model. effect of maternal micg vs current metsch. early life
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=20, alpha=0.05)
###restricted mean model. effect of maternal micg vs current metsch. total life
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=NULL, alpha=0.05)
#results indicate early life difference. no total life difference.

#recode tgvormet into new 0/1 variables with anone as 0. metsch as 1.
tgvvm2$arm<-tgvvm2$tgvormet
tgvvm2$arm<-recode_factor(tgvvm2$arm, anone="0")
tgvvm2$arm<-recode_factor(tgvvm2$arm, met="1")
###restricted mean model. effect of no exposure vs current metsch
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=20, alpha=0.05)
###restricted mean model. effect of no exposure vs current metsch
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=NULL, alpha=0.05)
#results indicate no early life difference. total life difference.

#recode tgvormet into new 0/1 variables with anone as 0. tgv as 1.
tgvvm2$arm<-tgvvm2$tgvormet
tgvvm2$arm<-recode_factor(tgvvm2$arm, anone="0")
tgvvm2$arm<-recode_factor(tgvvm2$arm, tgv="1")
###restricted mean model. effect of no exposure vs maternal micg. early life
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=20, alpha=0.05)
###restricted mean model. effect of no exposure vs maternal micg. total life
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=NULL, alpha=0.05)
#results indicate early life differnce. no total life difference.



lifesbadplot<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ tgvormet, data = tgvvm2),
                      xlab = "Days", 
                      ylab = "Survival Probability",
                      conf.int = FALSE, 
                      risk.table = FALSE,
                      censor=FALSE,
                      #linetype = c("solid","solid","dashed","dashed"),
                      size = 1, 
                      xlim=c(0,80),
                      ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                      #palette = c("BLACK", "GREY", "BLACK", "GREY"),
                      legend = c(.2,.1), 
                      legend.title = "Parasite Exposure")
                      #legend.labs = c( "None/None", "None/MicG","MicG/None", "MicG/MicG"))


lifesbadplot


#remove row 20 as this animal was not infected with metsch
# remove rows in r by row number
tgvvm2 <- tgvvm[-c(20),]

###restricted mean model. effect of maternal micg vs current metsch
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=20, alpha=0.05)
###restricted mean model. effect of maternal micg vs current metsch
rmst2(tgvvm2$deathday, tgvvm2$died, tgvvm2$arm, tau=78, alpha=0.05)

lifesbadplot2<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ tgvormet, data = tgvvm2),
                         xlab = "Days", 
                         ylab = "Survival Probability",
                         conf.int = FALSE, 
                         risk.table = FALSE,
                         censor=FALSE,
                         #linetype = c("solid","solid","dashed","dashed"),
                         size = 1, 
                         xlim=c(0,80),
                         ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                         #palette = c("BLACK", "GREY", "BLACK", "GREY"),
                         legend = c(.2,.1), 
                         legend.title = "Parasite Exposure")
#legend.labs = c( "None/None", "None/MicG","MicG/None", "MicG/MicG"))


lifesbadplot2

