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

redo$Parasite<-recode_factor(redo$Parasite, MicG="O. pajunii")
neosize$Parasite<-recode_factor(neosize$mpara, MicG="O. pajunii")

supfig$MomPXP<-paste(supfig$MomParasite,supfig$Parasite)

PanelA<-ggsurvplot(fit=survfit(Surv(deathday,died)~mpXpXc, data = supfig),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   censor=FALSE,
                   linetype = c("dotted", "solid", "dotted","solid"),
                   size = 2, 
                   facet_grid(.~Clone),
                   xlim=c(0,20),
                   ggtheme = theme_classic2(base_size=18, base_family = "Arial"),
                   palette = c("BLACK","BLACK" , "#898989", "#898989"),
                   legend = c(.25,.25), 
                   legend.title = "Maternal Exp/Current Exposure")
                   #legend.labs = c("O.p/O.p", "O.p/None", "None/O.p", "None/None"))+
                  #labs(tag="A")  #label for multipanel figure
    
  
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



PanelC<-ggplot(data=redo, aes(x = Parasite, y=sum, fill = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values =c("#1F1F1F","#898989","")) +
  #scale_pattern_manual(values = c(Exposed = "circle", None = "none")) +
  labs(x = "", y = "Early Reproduction", fill = "Maternal Exposure") + 
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))  +  
labs(tag="C")#label for multipanel figure

PanelC

####NUMBER OF CLUTCHES

PanelD<-ggplot(data=redo, aes(x = Parasite, y=clutches, fill = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values =c("#1F1F1F","#898989","")) +
  #scale_pattern_manual(values = c(Exposed = "circle", None = "none")) +
  labs(x = "Current Parasite", y = "Clutches Produced", fill = "Maternal Exposure") + 
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17)) +   
labs(tag="D")#label for multipanel figure

PanelD


library(ggplot2)
library(dplyr)

PanelB<-ggplot(data=neosize, aes(x = Parasite, y=Size, fill = Gmapara)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values =c("#1F1F1F","#898989","")) +
  labs(x = "", y = "Neonate Size", fill = "Maternal Exposure") + 
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17)) +   
  labs(tag="B")#label for multipanel figure

PanelB

#Making the two panel plot 

Panela<-PanelA$plot
grid.arrange(arrangeGrob(Panela,ncol=1,nrow=1),arrangeGrob( PanelB, ncol=1,nrow=1),arrangeGrob(PanelC, ncol=1,nrow=1), arrangeGrob(PanelD, ncol=1,nrow=1),widths=c(2))

TGVf1<-grid.arrange(Panela, arrangeGrob(PanelB, PanelC, PanelD), ncol = 2, widths=c(2,1))
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

