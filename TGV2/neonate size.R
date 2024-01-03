###analyses of neonate size for s, b,d clonal genotypes
###data used: TGV neonate size 
library(dplyr)
library(lme4)
neonatesize$Para<-recode_factor(neonatesize$Para, None="aNone")
neonatesize$GMApara<-recode_factor(neonatesize$GMApara, None="aNone")
neonatesize$GMAxMpara<-paste(neonatesize$GMApara,neonatesize$Para)

#create dataframe for each clone
Sneonatesize<-neonatesize[neonatesize$Clone=="S",]
Dneonatesize<-neonatesize[neonatesize$Clone=="D",]
Bneonatesize<-neonatesize[neonatesize$Clone=="B",]

#create only current unexposed dataset
unSneonatesize<-neonatesize[neonatesize$Para=="aNone",]
unDneonatesize<-neonatesize[neonatesize$Para=="None",]
unBneonatesize<-neonatesize[neonatesize$Para=="None",]


#####neonate size S clone
#Sneoglm <- glm(gdaughtersize ~ GMApara*Para + (GMApara*Para/MomID),data=Sneonatesize)
Sneolm <- lmer(gdaughtersize ~ GMApara*Para + (1|MomID),data=Sneonatesize)
Sneolm <- lmer(gdaughtersize ~ GMApara * (1|MomID),data=unSneonatesize)
summary(Sneolm)
anova(Sneolm, test="F")
###check the random effects
ranef(Sneolm)
#test linearity assumption. looking for randomly assorted values. plots resids, observed vs. predicted
Plot.Sneoglm.Linearity<-plot(resid(Sneolm),Sneonatesize$gdaughtersize)
#test homogeneity of variances. p>0.05 indicates no differences in variances
Sneonatesize$Model.S.Res<- residuals(Sneolm) #extracts the residuals and places them in a new column in our original data table
Sneonatesize$Abs.Model.S.Res <-abs(Sneonatesize$Model.S.Res) #creates a new column with the absolute value of the residuals
Sneonatesize$Model.S.Res2 <- Sneonatesize$Abs.Model.S.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.Model.S <- lm(Model.S.Res2 ~ MomID, data=Sneonatesize) #ANOVA of the squared residuals
anova(Levene.Model.S) #displays the results
##check normality
#id: identifies values that may be exerting undue influence on the model (i.e. outliers)
library(lattice)
qqmath(Sneolm, id=0.05) 
#not normal. opt for nonpara.

library(ARTool)
Sneonatesize$MomID<-as.factor(Sneonatesize$MomID)
Sneonatesize$GMApara<-as.factor(Sneonatesize$GMApara)
Sneonatesize$Para<-as.factor(Sneonatesize$Para)
Sneo<-art(gdaughtersize ~GMApara*Para, data = Sneonatesize)
summary(Sneo)
anova(Sneo)

art.con(Sneo, "GMApara:Para")

Sneonatesize$GmaXma<-paste(Sneonatesize$GMApara,Sneonatesize$Para)
aggregate(Sneonatesize$gdaughtersize, list(Sneonatesize$GmaXma), FUN=mean)
aggregate(Sneonatesize$gdaughtersize, list(Sneonatesize$GmaXma), FUN=sd)
aggregate(Sneonatesize$gdaughtersize, list(Sneonatesize$Para), FUN=mean)

#####neonate size B clone
Bneoglm <- lmer(gdaughtersize ~ GMApara*Para + (1|MomID),data=Bneonatesize)
summary(Bneoglm)
anova(Bneoglm, test="F")
##check residuals. 
##for variance= red lines should be centered on 0 or 1 and horiz.
##for normality qq should should align on diag.
par(mfrow=c(2,2))
plot(Bneoglm)
par(mfrow=c(1,1))
library(lattice)
qqmath(Bneoglm, id=0.05) 
#not normal. opt for nonpara.

library(ARTool)
Bneonatesize$MomID<-as.factor(Bneonatesize$MomID)
Bneonatesize$GMApara<-as.factor(Bneonatesize$GMApara)
Bneonatesize$Para<-as.factor(Bneonatesize$Para)
Bneo<-art(gdaughtersize ~GMApara*Para, data = Bneonatesize)
summary(Bneo)
anova(Bneo)

#art.con(Bneo, "GMApara:Para")
#Sneonatesize$GmaXma<-paste(Sneonatesize$GMApara,Sneonatesize$Para)
aggregate(Bneonatesize$gdaughtersize, list(Bneonatesize$GMApara), FUN=mean)
aggregate(Bneonatesize$gdaughtersize, list(Bneonatesize$GMApara), FUN=sd)
aggregate(Bneonatesize$gdaughtersize, list(Bneonatesize$Para), FUN=mean)
aggregate(Bneonatesize$gdaughtersize, list(Bneonatesize$Para), FUN=sd)



#####neonate size D clone
Dneoglm <- lmer(gdaughtersize ~ GMApara*Para + (1|MomID),data=Dneonatesize)
summary(Dneoglm)
anova(Dneoglm, test="F")
##check residuals. 
##for variance= red lines should be centered on 0 or 1 and horiz.
##for normality qq should should align on diag.
par(mfrow=c(2,2))
plot(Dneoglm)
par(mfrow=c(1,1))
##normality. 
library(lattice)
qqmath(Sneolm, id=0.05) 
#not normal. opt for nonpara.
###raw and transformed size does not fit well


library(ARTool)
Dneonatesize$MomID<-as.factor(Dneonatesize$MomID)
Dneonatesize$GMApara<-as.factor(Dneonatesize$GMApara)
Dneonatesize$Para<-as.factor(Dneonatesize$Para)
Dneo<-art(gdaughtersize ~GMApara*Para, data = Dneonatesize)
summary(Dneo)
anova(Dneo)

art.con(Dneo, "GMApara:Para")
Dneonatesize$GmaXma<-paste(Dneonatesize$GMApara,Dneonatesize$Para)
aggregate(Dneonatesize$gdaughtersize, list(Dneonatesize$GmaXma), FUN=mean)
aggregate(Dneonatesize$gdaughtersize, list(Dneonatesize$GmaXma), FUN=sd)

PanelD<-ggplot(data=neonatesize, aes(x = GMApara, y=gdaughtersize, fill = Para)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values = colorRampPalette(c("WHITE","BLACK"))(2)) +
  labs(x = "Maternal Exposure", y = "Neonate size") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="top")+ theme(text = element_text(size = 17)) +   
facet_grid(. ~ Clone)

PanelD


#####neonate size S clone average
Saveneoglm <- glm(neomean ~ GMApara*Para,data=Sneonatesize)
summary(Saveneoglm)
anova(Saveneoglm, test="F")
##check residuals. 
##for variance= red lines should be centered on 0 or 1 and horiz.
##for normality qq should should align on diag.
par(mfrow=c(2,2))
plot(Saveneoglm)
par(mfrow=c(1,1))
##normality. p<0.05=violation of assump.
shapiro.test(rstandard(Saveneoglm))
###raw and transformed size does not fit well

Saveneosize <- na.omit(Sneonatesize)
Saveneo<-art(neomean ~GMApara*Para, data = Saveneosize)
summary(Saveneo)
anova(Saveneo)
