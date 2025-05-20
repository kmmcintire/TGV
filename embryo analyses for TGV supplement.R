###EMBRYO EXPERIMENT ANALYSES
#Analyses included:
#Figure S1 created
#KMM 5.19.25

#Datafile in use: embryo_5.19.xlsx
#SHEETS:
#edata imported as edata
#neo_weight imported as wt
#ranges imported as ranges
#embryo_averages imported as aves



library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(forcats)
library(ggpubr)
library(car)
library(ARTool)
library(gridExtra)



###create subsets from edata based on treatment
#a=embryos for weighing, not allowed to develop
#b=embryos removed to develop in dish
#c=embryos developing in mother
justc<-edata[edata$TX=="Not Removed",]
justb<-edata[edata$TX=="Removed",]
justbc<-edata[edata$TX=="Removed"| edata$TX=="Not Removed",]

#mean and sd # of embryo by mom parasite
edata2<-edata[!is.na(edata$neonates),]
edata %>%
  group_by(momp) %>%
  get_summary_stats(embryos, type = "mean_sd")

#mean and sd # of neonates by mom parasite
edata %>%
  group_by(momp) %>%
  get_summary_stats(neonates, type = "mean_sd")

#mean and sd difference in embryo produced - neonate produced by mom parasite
edata %>%
  group_by(momp) %>%
  get_summary_stats(difference, type = "mean_sd")

#mean and sd ruptured embryos within the difference count by mom parasite
edata %>%
  group_by(momp) %>%
  get_summary_stats(difcountrupt, type = "mean_sd")


###in B treatment, does the number of neonates differ from embryos 
###the same for both groups
res <- t.test(difference ~ momp, data = justb)
res
###in B treatment, does number of neonates differ from embryos,
## the same for both groups if all rupture is handler error
res <- t.test(difcountrupt ~ momp, data = justb)
res

res <- t.test(embryos ~ momp, data = justb)
res

res <- t.test(neonates ~ momp, data = justb)
res
###THESE ANALYSES INDICATE:
# no difference in embryos produced but a difference in survival of those animals dependent on momparasite

###in C treatment, does the number of neonates differ from embryos 
###the same for both groups
res <- t.test(difference ~ momp, data = justc)
res
###in c treatment, does number of neonates differ from embryos,
## the same for both groups if all rupture is handler error
res <- t.test(difcountrupt ~ momp, data = justc)
res

res <- t.test(embryos ~ momp, data = justc)
res

res <- t.test(neonates ~ momp, data = justc)
res
####no differences dependent on momparasite


###Analysis of neonate #
res.aov2 <- aov(neonates ~ momp * TX, data = justbc)
summary(res.aov2)
#assess homogeneity ofvariances
plot(res.aov2, 1)
leveneTest(neonates ~ momp * TX , data = justbc)
#assess Normality assumption
plot(res.aov2, 2)
# Extract the residuals
aov_residuals <- residuals(object = res.aov2)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )
#data do not meet assumptions. opt for nonparametric analysis

####USED IN MANUSCRIPT, SUPPLEMENT. Analysis for neonates produced. Maternal care TX B and C.
justbc$momp<-as.factor(justbc$momp)
justbc$TX<-as.factor(justbc$TX)
justbc <- justbc %>% drop_na(embryos)  
neo = art(neonates ~ momp * TX , data=justbc) # linear mixed model syntax; see lme4::lmer
summary(neo)
anova(neo)
art.con(neo, "momp:TX")
##Indicates significant interaction

###PANEL B FOR SUPPLEMENTAL FIGURE 1#########################
justbc$Matpara<-recode_factor(justbc$Matpara, control = "Control")
panelB <- justbc %>%
  mutate(
    Matpara = fct_relevel(Matpara, "O.pajunii", "Control"),
    EmbryoCare = fct_relevel(EmbryoCare, "Removed", "Not Removed")  
  ) %>%
  ggplot( aes(x = Matpara, y = neonates, color=EmbryoCare)) + 
  geom_boxplot(lwd=1)+
  labs(title="",x="", y = "Neonate Abundance")+
  #scale_fill_manual(values=c("#898989", "white"))+
  scale_color_manual(values=c( "#898989","black"))+
  theme_classic()+
  theme(axis.title = element_text(face="bold", colour="black", size=12))+
  theme(legend.position="top")
panelB





####USED IN MANUSCRIPT, supplement. ANALYSIS of difference between neonates and embryos (embryo loss), includes b and c
dif = art(difference ~ momp * TX , data=justbc) 
summary(dif)
anova(dif)
art.con(dif, "momp:TX")
##significant interaction of effects

####panel A for manu figure, supplement figure 1############
panelA <- justbc %>%
  mutate(
    Matpara = fct_relevel(Matpara, "O.pajunii", "Control"),
    EmbryoCare = fct_relevel(EmbryoCare, "Removed", "Not Removed")  
  ) %>%
  ggplot(aes(x = Matpara, y = difference, color=EmbryoCare)) + 
  geom_boxplot(lwd=1)+
  labs(title="",x="", y = "Embryos Lost")+
  #scale_fill_manual(values=c("#898989", "white"))+
  scale_color_manual(values=c("#898989", "black"))+
  theme_classic()+
  theme(axis.title = element_text(face="bold", colour="black", size=12))+
  theme(legend.position="")
panelA


#####USED IN MANUSCRIPT, supplement. embryo production.
#no need for embryo care tx, as this had yet not occurred at time of embryo production
res <- t.test(embryos ~ momp, data = edata)
res
# Shapiro-Wilk normality test for control's sizes
with(edata, shapiro.test(embryos[momp == "control"]))
# Shapiro-Wilk normality test for micg's sizes
with(edata, shapiro.test(embryos[momp == "MicG"]))
##nonnormality. opt for wilcoxon
###micg not normal, opt for wilcoxon
res <- wilcox.test(embryos ~ momp, data = edata,
                   exact = FALSE)
res
###No difference in embryo produced for all treatments
ggboxplot(edata, x = "momp", y = "embryos", color = "momp",
          palette = c("#00AFBB", "#E7B800") )


###USED IN MANUSCRIPT supplement. neonate average mass (indave). tx C only.
res <- t.test(indave ~ matpara, data = wt)
res
# Shapiro-Wilk normality test for control's sizes
with(wt, shapiro.test(indave[matpara == "C"]))
# Shapiro-Wilk normality test for micg's sizes
with(wt, shapiro.test(indave[matpara == "O.pajunii"]))
##normality confirmed.
###no difference based on maternal pathogen



####USED IN MANUSCRIPT supplement. is tehre greater variation within mom dependent on mom pathogen exposure?
##embryo size
#ONLY embryo care tx A
res <- t.test(sizerange ~ momp, data = ranges)
res
# Shapiro-Wilk normality test for control's sizes
with(ranges, shapiro.test(sizerange[momp == "Control"]))
# Shapiro-Wilk normality test for micg's sizes
with(ranges, shapiro.test(sizerange[momp == "O.pajunii"]))
##nonnormality. opt for wilcoxon
###micg not normal, opt for wilcoxon
res <- wilcox.test(sizerange ~ momp, data = ranges,
                   exact = FALSE)
res

###FIG for manu, SUpplemental FIgure 1 panel D###############

panelD<-ranges %>%
  mutate(
    momp = fct_relevel(momp, "O.pajunii", "Control")
    )%>%  
  ggplot( aes(x = momp, y = sizerange)) + 
  geom_boxplot(lwd=1)+
  labs(title="",x="Maternal Exposure", y = "Within Clutch Embryo \n Length Range")+
  #scale_fill_manual(values=c("#898989", "white"))+
  scale_color_manual(values=c("black", "black"))+
  theme_classic()+
  theme(axis.title = element_text(face="bold", colour="black", size=12))+
  theme(legend.position="")
panelD





######USED IN MANUSCRIPT supplement. is there greater within clutch variation in yolk size dependent on momparasite exposure?
#yolk size range
###ONLY embryo care tx A
res <- t.test(yolkrange ~ momp, data = ranges)
res
# Shapiro-Wilk normality test for control's sizes
with(ranges, shapiro.test(yolkrange[momp == "Control"]))
# Shapiro-Wilk normality test for micg's sizes
with(ranges, shapiro.test(yolkrange[momp == "O.pajunii"]))
##nonnormality. opt for wilcoxon
###micg not normal, opt for wilcoxon
res <- wilcox.test(yolkrange ~ momp, data = ranges,
                   exact = FALSE)
res


#####FIG for manu Supplemental FIgure 1 panel C #####################
panelC<-ranges %>%
  mutate(
    momp = fct_relevel(momp, "O.pajunii", "Control")
  )%>%  
  ggplot(aes(x = momp, y = yolkrange )) + 
  geom_boxplot(lwd=1)+
  labs(title="",x="Maternal Exposure", y = "Within Clutch \n Yolk Width Range")+
  #scale_fill_manual(values=c("#898989", "white"))+
  scale_color_manual(values=c("black", "black"))+
  theme_classic()+
  theme(axis.title = element_text(face="bold", colour="black", size=12))+
  theme(legend.position="")
panelC



###########USED IN MANUSCRIPT supplement. Analysis of embryonic length, embryos from each mom averaged
res <- t.test(size ~ momp, data = aves)
res
# Shapiro-Wilk normality test for control's sizes
with(embryo, shapiro.test(size[momp == "Control"]))
# Shapiro-Wilk normality test for micg's sizes
with(embryo, shapiro.test(size[momp == "MicG"]))
#test for homogen variances.
res.ftest <- var.test(size ~ momp, data = embryo)
res.ftest #does not pass, not required for welchs ttest

########USED IN MANUSCRIPT supplement.Analysis of yolk/length ratio, avereaged by mom
res <- t.test(yolkratio ~ momp, data = aves)
res
# Shapiro-Wilk normality test for control's sizes
with(embryo, shapiro.test(yolkratio[momp == "Control"]))
# Shapiro-Wilk normality test for micg's sizes
with(embryo, shapiro.test(yolkratio[momp == "MicG"]))

#######USED IN MANUSCRIPT supplement. analysis of yolk width, averaged by mom
res <- t.test(yolk ~ momp, data = aves)
res
# Shapiro-Wilk normality test for control's sizes
with(embryo, shapiro.test(yolk[momp == "Control"]))
# Shapiro-Wilk normality test for micg's sizes
with(embryo, shapiro.test(yolk[momp == "MicG"]))
###micg yolk sizes are not normally dist. opt for wilcox
res.ftest <- var.test(yolk ~ momp, data = embryo)
res.ftest
res <- wilcox.test(yolk ~ momp, data = aves,
                   exact = FALSE)
res



######USED IN MANUSCRIPT supplement.Anaysis of embryosize/momsize ratio. averaged by mom
res <- t.test(momratio ~ momp, data = aves)
res
# Shapiro-Wilk normality test for control's sizes
with(embryo, shapiro.test(momratio[momp == "Control"]))
# Shapiro-Wilk normality test for micg's sizes
with(embryo, shapiro.test(momratio[momp == "MicG"]))
##nonnormality. opt for wilcoxon
res <- wilcox.test(momratio ~ momp, data = aves,
                   exact = FALSE)
res



####PUT TOGETHER SUPPLEMENTAL FIGURE 1
grid.arrange(arrangeGrob(panelA,ncol=1,nrow=1),arrangeGrob(panelB, ncol=1,nrow=1),arrangeGrob(panelC, ncol=1,nrow=1), widths=c(2) )

TGV_SUPFIG1<-ggarrange(panelA, panelB, panelC,panelD, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)


TGV_SUPFIG1
ggsave("TGV_SUPFIG1.png", TGV_SUPFIG1)


