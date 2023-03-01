
###Making the TGV survival figure
install.packages("extrafont")
library(extrafont)
font_import()
loadfonts()

install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("cli")
install.packages("tidyverse")
library(tidyverse)
library(cli)
library(ggplot2)
library(survival)
library(ggpubr)
library(dplyr)
library(survminer)
install.packages("gridExtra")
library(gridExtra)

controlbabies$Clone<-recode_factor(controlbabies$Clone, STD="S")
controlbabies$Clone<-recode_factor(controlbabies$Clone, BD="B")
controlbabies$Clone<-recode_factor(controlbabies$Clone, DW="D")

PanelA<-ggsurvplot(fit=survfit(Surv(deathday,censor)~MomPXClone, data = controlbabies),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   censor=FALSE,
                   linetype = c("dotted", "dotted", "dotted","solid","solid","solid"),
                   size = 2, 
                   xlim=c(0,20),
                   ggtheme = theme_classic2(base_size=18, base_family = "Arial"),
                   palette = c("#5e3c99", 
                               "#1B9E77", "#e66101", "#5e3c99", 
                               "#1B9E77", "#e66101"), 
                   legend = c(.2,.25), 
                   legend.title = "",
                   legend.labs = c("M B", 
                                  "M D","M S", "N B", "N D", "N S"))


PanelA
####MAKING TGV REPRODUCTION FIGURE
install.packages("MASS")
library(MASS)
install.packages("sjPlot")
library(sjPlot)
install.packages("sjmisc")
library(sjmisc)
install.packages("ggplot2")
library(ggplot2)
install.packages("ggpattern")
library(ggpattern)
library(dplyr)

controlbabies$Clone<-as.character(controlbabies$Clone)

PanelB<-ggplot(data=controlbabies, aes(x = Clone, y=sumrepro, fill = Clone, pattern = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  geom_boxplot_pattern(position = position_dodge(preserve = "single"),
                                                        color = "black", 
                                                        pattern_fill = "black",
                                                        pattern_angle = 45,
                                                        pattern_density = 0.5,
                                                        pattern_spacing = 0.025,
                                                        pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = colorRampPalette(c("#5e3c99","#1B9E77","#e66101"))(3)) +
  scale_pattern_manual(values = c(Exposed = "circle", None = "none")) +
  labs(x = "Genotype", y = "Total Reproduction", pattern = "Maternal Exposure") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))    
                                                    
                                                                                                            
PanelB
library(ggplot2)
library(dplyr)
babysize$Clone<-recode_factor(babysize$Clone, STD="S")
babysize$Clone<-recode_factor(babysize$Clone, DW="D")
babysize$Clone<-recode_factor(babysize$Clone, BD="B")

PanelC<-ggplot(data=babysize, aes(x = Clone, y=Size, fill = Clone, pattern = GMApara)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  geom_boxplot_pattern(position = position_dodge(preserve = "single"),
                                                        color = "black", 
                                                        pattern_fill = "black",
                                                        pattern_angle = 45,
                                                        pattern_density = 0.5,
                                                        pattern_spacing = 0.025,
                                                        pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = colorRampPalette(c("#5e3c99","#1B9E77","#e66101"))(3)) +
  scale_pattern_manual(values = c(MicG = "circle", None = "none")) +
  labs(x = "Genotype", y = "Neonate Size", pattern = "Grandmaternal Exposure") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")   + theme(text = element_text(size = 17)) 
PanelC


#Making the two panel plot 

Panela<-PanelA$plot
grid.arrange(arrangeGrob(Panela,ncol=1,nrow=1),arrangeGrob( PanelB, ncol=1,nrow=1),arrangeGrob(PanelC, ncol=1,nrow=1), widths=c(2) )

TGVf1<-grid.arrange(Panela, arrangeGrob(PanelB, PanelC), ncol = 2, widths=c(2,1))
ggsave("TGVfig1.png", TGVf1)

Panela<-Panelsoliddot$plot




#####Making the two legends needed to piece it together

###making the solid/dot legend
Panelsoliddot<-ggsurvplot(fit=survfit(Surv(deathday,censor)~MomParasite, data = controlbabies),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   linetype = c("dotted","solid"),
                   size = 1, 
                   censor= FALSE,
                   xlim=c(0,20),
                   ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                   palette = c("BLACK", 
                               "BLACK"), 
                   legend = c(.25,.25),
                   legend.title = "Maternal Parasite",
                   legend.labs = c("Exposed", 
                                   "None"))


Panelsoliddot

PanelA<-ggsurvplot(fit=survfit(Surv(deathday,censor)~Clone, data = controlbabies),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   linetype = c("solid","solid","solid"),
                   size = 1, 
                   censor= FALSE,
                   xlim=c(0,20),
                   ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                   palette = c("#5e3c99", 
                               "#1B9E77", "#e66101"), 
                   legend = c(.25,.25),
                   legend.title = "Genotype",
                   legend.labs = c("B", 
                                   "D","S"))


PanelA

PanelA<-ggsurvplot(fit=survfit(Surv(deathday,censor)~MomPXClone, data = controlbabies),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   linetype = c("solid","solid","solid", "dotted", "dotted", "dotted"),
                   size = 1, 
                   xlim=c(0,20),
                   ggtheme = theme_classic2(base_size=20, base_family = "Arial"),
                   palette = c("#5e3c99", 
                               "#1B9E77", "#e66101", "#5e3c99", 
                               "#1B9E77", "#e66101"), 
                   legend = c("none"))


PanelA

