####TGV SUPPLEMENTAL FIGURE WITH 9 panels
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

###bottom row: NEONATE SIZE
##data: neonatesize. Clone, GMApara, Para
neonatesize$Para<-recode_factor(neonatesize$Para, aNone="None")
neonatesize$GMApara<-recode_factor(neonatesize$GMApara, aNone="None")
library(dplyr)
library(forcats)
bottomrow<-ggplot(data=neonatesize, aes(x = GMApara, y=gdaughtersize, fill = Para)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values =c("#1F1F1F","#898989","")) +
  #scale_pattern_manual(values = c(Exposed = "circle", None = "none")) +
  labs(x = "Maternal Exposure", y = "Neonate Size", fill = "Exposure") + 
  facet_grid(.~Clone)+ #separates into facets by clone
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))  +  
  labs(tag="C")#label for multipanel figure

bottomrow
###

###Top ROW: c1-3 offspring #
###data tgvtryagain, sheet2. MomParasite, Parasite,Clone, sum


toprow<-ggplot(data=sheet2, aes(x = MomParasite, y=sum, fill = Parasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values =c("#1F1F1F","#898989","")) +
  #scale_pattern_manual(values = c(Exposed = "circle", None = "none")) +
  labs(x = "", y = "Early Reproduction", fill = "Exposure") + 
  facet_grid(.~Clone)+ #separates into facets by clone
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+  theme(text = element_text(size = 17))  +  
  labs(tag="A")#label for multipanel figure

toprow=toprow+ theme(legend.position=c(x=.5, y=.8), legend.background=element_blank(), legend.title=element_blank(), legend.box.background = element_rect(colour="black") ,legend.text = element_text(size = 6), legend.direction="horizontal")
toprow

###Middle row: # of clutches
###data tgvtryagain, sheet2. MomParasite, Parasite, clutches

middlerow<-ggplot(data=sheet2, aes(x = MomParasite, y=clutches, fill = Parasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot()+theme_classic()+  
  scale_fill_manual(values =c("#1F1F1F","#898989","")) +
  #scale_pattern_manual(values = c(Exposed = "circle", None = "none")) +
  labs(x = "", y = "Clutches Produced",fill = "Exposure") + 
  facet_grid(.~Clone)+ #separates into facets by clone
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))  +  
  labs(tag="B")#label for multipanel figure

middlerow
sfig<-grid.arrange(toprow,middlerow, bottomrow, ncol = 1)
sfig

