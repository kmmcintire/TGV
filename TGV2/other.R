
supfig$Parasite<-recode_factor(supfig$Parasite, None="aNone")
supfig$MomParasite<-recode_factor(supfig$MomParasite, None="aNone")
supfig$MomPXClone<-paste(supfig$MomParasite,supfig$Clone)
supfig$mpXpXc<-paste(supfig$MomPXClone, supfig$Parasite)

supfig$Clone<-recode_factor(supfig$Clone, BD="B")
supfig$Clone<-recode_factor(supfig$Clone, DW="D")
supfig$Clone<-recode_factor(supfig$Clone, STD="S")

lifeall<-coxph(Surv(deathday,edied)~ MomParasite + Parasite * Clone, data=supfig)
summary(lifeall)
test.ph<-cox.zph(lifeall)
test.ph
par(mfrow=c(2,2))
plot(cox.zph(lifeall))

PanelA<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomPXPara, data = redoB),
                      xlab = "Days", 
                      ylab = "Survival Probability",
                      conf.int = FALSE, 
                      risk.table = FALSE,
                      censor=FALSE,
                      linetype = c("solid","dashed","solid","dashed"),
                      size = 1.5, 
                      xlim=c(0,21),
                      break.x.by = 10,
                      #show.legend=FALSE,
                      ggtheme = theme_classic2(base_size=14, base_family = "Arial"),
                      palette = c("#898989","#898989","BLACK","BLACK"),
                      legend = c(.25,.3),
                      legend.title = "Maternal/Current",
                   legend.text = element_text(size=12),
                      legend.labs = c( "None/None","None/O.p", "O.p/None", "O.p/O.p"))+
                        labs(tag="A")  #label for multipanel figure

#PanelA
PanelA$plot <- PanelA$plot+ 
  ggplot2::annotate("text", 
                    x = 16, y = 0.1, # x and y coordinates of the text
                    label = "B genotype \n RMST diff = -1.844 \n p = 0.083", size = 4)
PanelA

PanelB<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomPXPara, data = redoD),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   censor=FALSE,
                   linetype = c("solid","dashed","solid","dashed"),
                   size = 1.5, 
                   xlim=c(0,21),
                   break.x.by = 10,
                   show.legend=FALSE,
                   ggtheme = theme_classic2(base_size=14, base_family = "Arial"),
                   palette = c("#898989","#898989","BLACK","BLACK"),
                   legend = c(.25,.2),
                   legend.title = "Maternal/Current",
                   legend.labs = c( "None/None","None/O.pajunii", "O.pajunii/None", "O.pajunii/O.pajunii"))+
  labs(tag="B")  #label for multipanel figure

#PanelB
PanelB$plot <- PanelB$plot+ 
  ggplot2::annotate("text", 
                    x = 15, y = 0.1, # x and y coordinates of the text
                    label = "D genotype \n RMST diff = -4.3 \n p = 0.008", size = 4)
PanelB


PanelC<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomPXPara, data = redo),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   censor=FALSE,
                   linetype = c("solid","dashed","solid","dashed"),
                   size = 1.5, 
                   xlim=c(0,21),
                   break.x.by = 10,
                   show.legend=FALSE,
                   ggtheme = theme_classic2(base_size=14, base_family = "Arial"),
                   palette = c("#898989","#898989","BLACK","BLACK"),
                   legend = c(.25,.2),
                   legend.title = "Maternal/Current",
                   legend.labs = c( "None/None","None/O.p", "O.p/None", "O.p/O.p"))+
  labs(tag="C")  #label for multipanel figure

#PanelC
PanelC$plot <- PanelC$plot+ 
  ggplot2::annotate("text", 
                    x = 7, y = 0.1, # x and y coordinates of the text
                    label = "S genotype \n RMST diff = -4.141 \n p = 0.004", size = 4)
PanelC

PanelD<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomParasite, data = immact),
                   xlab = "Days", 
                   ylab = "Survival Probability",
                   conf.int = FALSE, 
                   risk.table = FALSE,
                   censor=FALSE,
                   #linetype = c("solid","dashed","solid","dashed","solid", "dashed","solid","dashed", "solid","dashed", "solid","dashed"),
                   size = 1.5, 
                   xlim=c(0,80),
                   break.x.by = 10,
                   ggtheme = theme_classic2(base_size=14, base_family = "Arial"),
                   palette = c("BLACK","#898989","BLUE","RED"),
                    legend = c(.90,.77),
legend.title = "Maternal",
#legend.title = element_text(size=10), #change legend title font size
legend.text = element_text(size=12), #change legend text font size
legend.labs = c( "O.p","None","P.r", "M.b"))+
labs(tag="D")  #label for multipanel figure

PanelD




###this panel creates visual for the TGVvMETSCH analysis
###THIS figure INCLUDES animal 508
PanelE<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ tgvormet, data = tgvvm2),
                         xlab = "Days", 
                         ylab = "Survival Probability",
                         conf.int = FALSE, 
                         risk.table = FALSE,
                         censor=FALSE,
                         linetype = c("solid","dashed","solid"),
                         size = 1.5, 
                         xlim=c(0,80),
                         break.x.by = 10,
                         ggtheme = theme_classic2(base_size=14, base_family = "Arial"),
                         palette = c("#898989", "RED", "BLACK"),
                         legend = c(.85,.75), 
                   #legend.title = element_text(size=10), #change legend title font size
                   legend.text = element_text(size=12), #change legend text font size
                         legend.title = ("Maternal/Current"),
                        legend.labs = c( "None/None", "None/M.b", "O.p/None"))+
                        labs(tag="E")  #label for multipanel figure

PanelE

#Making the two panel plot 



#######################################
library(gridExtra)
Paneld<-PanelD$plot
Panele<-PanelE$plot
rowone<-grid.arrange(PanelA$plot,PanelB$plot,PanelC$plot,nrow=1)
rowtwo<-grid.arrange(Paneld, Panele, nrow = 1)
ggarrange(rowone,rowtwo, ncol=1, nrow=2)

             



fig1lifeplot<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomPXClone, data = supfig),
                        xlab = "Days", 
                        ylab = "Survival Probability",
                        conf.int = FALSE, 
                        risk.table = FALSE,
                        censor=FALSE,
                        linetype = c("solid","solid","solid", "dashed","dashed", "dashed"),
                        size = 1, 
                        xlim=c(0,21),
                        break.x.by = 10,
                        facet.by = "Clone",
                        ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                        palette = c("BLACK", "BLACK", "BLACK","BLACK","BLACK","BLACK"))
                        #legend = c(.25,.2), 
                        #legend.title = "Maternal Parasite",
                        #legend.labs = c( "None B","NoneD", "None S", "O.pajunii B","O.pajuniiD","O.pajunii s"))


fig1lifeplot

fig1lifeplotlegend<-ggsurvplot(fit=survfit(Surv(deathday, died) ~ MomParasite, data = supfig),
                         xlab = "Days", 
                         ylab = "Survival Probability",
                         conf.int = FALSE, 
                         risk.table = FALSE,
                         censor=FALSE,
                         linetype = c("solid","dashed"),
                         size = 1, 
                         xlim=c(0,21),
                         break.x.by = 10,
                         #facet.by = "Clone",
                         ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                         palette = c("BLACK", "BLACK"),
                         legend = c(.25,.2), 
                         legend.title = "Maternal Parasite",
                         legend.labs = c( "None","O.pajunii"))


fig1lifeplotlegend
