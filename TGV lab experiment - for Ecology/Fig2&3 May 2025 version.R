#FIGURE 2 & 3 used in TGV manuscript. 8.2.24

#title: Transgenerational pathogen effects: Maternal pathogen exposure reduces offspring fitness
#Submitted to Ecology for review
#Neosize sheet imported as "NeoSizeA" --used in F3 panel D
#TGV vs Metsch sheet imported as "tgvvmz" --Used in Figure 2
#Main EXP sheet imported as "redoz" --used in Figure 2, FIgure 3 panels A-C

library(survival)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(survminer)
library(gridExtra)
library(ggpubr)


library(extrafont)
font_import()
loadfonts(device = "win")      
windowsFonts()   





####################
#Figure 2
################

## main life plot
#gather and configure data for plot
#tgvvmz dataset was used for tgv vs metsch comparison, the metsch animals will be used in this plot
#remove uninfected metsch animal
tgvvmz2 <- tgvvmz
# remove rows in r by row number # remove uninfected metsch animal
tgvvmz2 <- tgvvmz2[-c(20),]
metforfig2z <- subset(tgvvmz2, (tgvormet == 'met'))
metforfig2z<- metforfig2z[ -c(2,4,6:15, 18:19) ]
names(metforfig2z)[2] <- "MomParasite"
names(metforfig2z)[3] <- "Parasite"
#redoz dataset was used for TGV comparisons, these animals will be used in this plot
tgvforfig2z<-redoz
tgvforfig2z<- tgvforfig2z[ -c(2,4,6:17,20:22) ]
tgvforfig2z$Parasite<-recode_factor(tgvforfig2z$Parasite, aNone="None")
tgvforfig2z$Parasite<-recode_factor(tgvforfig2z$Parasite, MicG="O.pajunii")
tgvforfig2z$MomParasite<-recode_factor(tgvforfig2z$MomParasite, aNone="None")
tgvforfig2z$MomParasite<-recode_factor(tgvforfig2z$MomParasite, MicG="O.pajunii")
#combining sets and creating combined variiable for figure 2
lifefigsetz<-rbind(metforfig2z,tgvforfig2z)
lifefigsetz$MomPXPara<-paste(lifefigsetz$MomParasite,lifefigsetz$Parasite)


Figure2asof71824<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomPXPara, data =lifefigsetz),
                             xlab = "Days", 
                             ylab = "Survival Probability",
                             conf.int = FALSE, 
                             risk.table = FALSE,
                             censor=FALSE,
                             linetype = c("solid","solid","dashed","solid","dashed"),
                             size = 1.0, 
                             xlim=c(0,80),
                             palette = c("RED", "BLACK", "BLACK", "BLUE", "BLUE"),
                             legend = c(.85,.8), 
                             legend.title = "Pathogen Exposure",
                             legend.labs = c( "None/Mb", "None/None", "None/Op","Op/None", "Op/Op"))
Figure2asof71824$plot$layers <- c(annotate(geom = "rect", xmin = 0, xmax = 21, ymin = 0, ymax = 1, fill = "grey95"), Figure2asof71824$plot$layers)
ggpar(Figure2asof71824, 
      font.main = c(20, "bold"),
      font.x = c(20, "bold"),
      font.y = c(20, "bold"),
      font.caption = c(20, "bold"), 
      font.legend = c(16, "bold"), 
      font.tickslab = c(20, "bold"))


Figure2asof71824


#####################
#Figure 3
#####################

#3 panel figure about reproduction
tgvforfig3<-redoz[-c(2,4,6:7,12,14:22)]
tgvforfig3$Parasite<-recode_factor(tgvforfig3$Parasite, aNone="None")
tgvforfig3$Parasite<-recode_factor(tgvforfig3$Parasite, "MicG"="Op")
tgvforfig3$MomParasite<-recode_factor(tgvforfig3$MomParasite, aNone="None")
tgvforfig3$MomParasite<-recode_factor(tgvforfig3$MomParasite, MicG ="Op")


# panel a = CLUTCHES PRODUCED 

PA<-ggplot(data=tgvforfig3, aes(x = Parasite, y=clutches, fill = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot(color="GREY")+theme_classic(base_size = 20)+  
  scale_fill_manual(values = colorRampPalette(c("BLUE","BLACK"))(2)) +
  geom_point(color= "GREY", position=position_jitterdodge(seed=123))+ #seed makes jitter reproducible (added 5.9.25KMM)
  labs(x = "", y = "Clutches Produced") + 
  #guides(pattern = guide_legend(override.aes = list(fill = "white")),
        # fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="top")+ theme(text = element_text(size = 17))+  
    labs(tag="B")#label for multipanel figure
PA<-PA+guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))    
#PB=PB+ theme(legend.position=c(x=.5, y=.8), legend.background=element_blank(), legend.title=element_blank(), legend.box.background = element_rect(colour="black") ,legend.text = element_text(size = 6), legend.direction="horizontal")
PA




PB<-ggplot(data=tgvforfig3, aes(x = Parasite, y=sum, fill = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot(color="GREY")+theme_classic(base_size = 20)+  
  scale_fill_manual(values =c("BLUE","BLACK","")) +
  geom_point(color= "GREY", position=position_jitterdodge(seed=123))+#seed makes jitter reproducible (added 5.9.25KMM)
  labs(x = "", y = "Early Reproduction", fill = "Maternal Exposure") + 
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+  theme(text = element_text(size = 17))  +  
  labs(tag="C")#label for multipanel figure
PB<-PB+guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))    

PB

##panel c - neonate size
sizeforpostera<-NeosizeA
names(sizeforpostera)[4] <- "Para"
names(sizeforpostera)[3] <- "GMApara"
names(sizeforpostera)[5] <- "gdaughtersize"
sizeforpostera$Para<-recode_factor(sizeforpostera$Para, "O.pajunii"="Op")
sizeforpostera$GMApara<-recode_factor(sizeforpostera$GMApara, "O.pajunii"="Op")

PC<-ggplot(data=sizeforpostera, aes(x = Para, y=gdaughtersize, fill = GMApara)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot(outlier.shape=NA, color= "GREY")+theme_classic( base_size=20)+  
  geom_point(color= "GREY", position=position_jitterdodge(seed=123))+#seed makes jitter reproducible (added 5.9.25KMM)
  scale_fill_manual(values =c("BLUE","BLACK","")) +
  #scale_pattern_manual(values = c(Exposed = "circle", None = "none")) +
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+  theme(text = element_text(size = 17))  +
  labs(x = "", y = "Neonate Size", fill = "Maternal Exposure") +
  labs(tag="D")
#PC<-PC+guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))    
PC=PC+ theme(legend.position=c(x=.8, y=.9),legend.text = element_text(size = 10), legend.title = element_text(size = 10),legend.background=element_blank() )


PC

PD<-ggplot(data=tgvforfig3, aes(x = Parasite, y=c1size, fill = MomParasite)) +theme(text = element_text(family = "Times New Roman"))+
  geom_boxplot(color="GREY")+theme_classic(base_size = 20)+  
  scale_fill_manual(values =c("BLUE","BLACK","")) +
  geom_point(color= "GREY", position=position_jitterdodge(seed=123))+#seed makes jitter reproducible (added 5.9.25KMM)
  labs(x = "", y = "Size at 1st Reproduction", fill = "Maternal Exposure") + 
  guides(fill = guide_legend(override.aes = list(pattern = "none")))+  theme(text = element_text(size = 17))  +  
  labs(tag="A")#label for multipanel figure
PD<-PD+guides(fill = guide_legend(override.aes = list(pattern = "none")))+ theme(legend.position="none")+ theme(text = element_text(size = 17))    

PD

fig3asof71824<-grid.arrange(PD,PA,PB, PC, nrow = 2)

# Annotate the figure by adding a common labels
annotate_figure(fig3asof71824,
                bottom = text_grob("Current Exposure", color = "black",
                               size = 17)
                             )


fig3asof71824






