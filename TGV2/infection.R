



micgbabies$INFtot<-micgbabies$INF
micgbabies$INFtot[is.na(micgbabies$INFtot)] <- 0


####infection yes/no:
####comparing models with and without interaction. not better with interaction.
inf1 <- glm(INF ~ MomParasite * Clone, family = binomial, data = micgbabies)
summary(inf1)
inf2 <- glm(INF ~ MomParasite + Clone, family = binomial, data = micgbabies)
summary(inf2)
anova(inf1, inf2, test = "Chisq")
anova(inf1, test = "Chisq")


####infection yes/no:
####comparing models with and without interaction. not better with interaction.
###inf2tot used in tgv manuscript
inf1tot <- glm(INFtot~ MomParasite * Clone, family = binomial, data = micgbabies)
summary(inf1tot)
inf2tot <- glm(INFtot ~ MomParasite + Clone, family = binomial, data = micgbabies)
summary(inf2tot)
anova(inf1tot, inf2tot, test = "Chi")
anova(inf1tot, test = "Rao")



i<-ggplot(dft1b, aes( y=Clone,fill=MomParasite , pattern = Clone)) +
  geom_bar()+theme_classic()+  geom_bar_pattern(position = position_dodge(preserve = "single"),
                                                        color = "black", 
                                                        pattern_fill = "black",
                                                        pattern_angle = 45,
                                                        pattern_density = 0.5,
                                                        pattern_spacing = 0.025,
                                                        pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = c(BD ="#5e3c99",DW="#b2abd2",STD="#e66101")) +
  scale_pattern_manual(values = c(MicG = "stripe", None = "none")) +
  labs(x = "Clone", y = "Frequency of Infection",fill="Clone", pattern = "Maternal Exposure") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
i



t1<- table(micgbabies$Clone, micgbabies$MomParasite, paste(micgbabies$Clone, micgbabies$MomParasite), 
                micgbabies$INFtot, dnn = c("Clone","MomParasite","Clone+Parasite","Infection"))
t1
ftable(t1)

t2<- table(micgbabies$Clone, micgbabies$MomParasite, paste(micgbabies$Clone, micgbabies$MomParasite), 
           micgbabies$INF, dnn = c("Clone","MomParasite","Clone+Parasite","Infection"))
t2
ftable(t1)

dft1<-as.data.frame.table(t1)
dft1b <- t1 %>%        
 as_tibble() %>%             
 tidyr::uncount(n) %>%              
 mutate_all(as.factor)


ggplot(data = dft1b) +
 geom_bar(mapping = aes(x = Infection))

ggplot(data = dft1b, mapping = aes(x=`Clone+Parasite`, fill=Infection))+
 geom_bar(position="dodge") +  
   labs(title = "Hair Color", 
        subtitle = "592 Statistics Students",
        caption = "(From R's built in HairEyeColor sample dataset)",
        y = "Number of Students", x = NULL)  

ggplot(data = dft1b) + 
  geom_bar(mapping = aes(`Clone+Parasite`, fill = Infection))


dft2b <- t2 %>%        
  as_tibble() %>%             
  tidyr::uncount(n) %>%              
  mutate_all(as.factor)


ggplot(data = dft2b, aes(`Clone+Parasite`, fill = Clone, pattern=Infection, color=Clone) ) + 
  theme_classic() +
  geom_bar() +
  geom_bar_pattern(color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = c(BD ="#5e3c99",DW="#b2abd2",STD="#e66101")) +
  scale_pattern_manual(values = c("1" = "stripe", "0" = "none")) +
  labs(x = "Clone + Maternal Exposure", y = "Infected Individuals",fill="Clone", pattern = "Infection") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))


ggplot(data = dft1b, aes(`Clone+Parasite`, fill = Clone, pattern=Infection, color=Clone) ) + 
  theme_classic() +
  geom_bar() +
  geom_bar_pattern(color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.5,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = c(BD ="#5e3c99",DW="#b2abd2",STD="#e66101")) +
  scale_pattern_manual(values = c("1" = "stripe", "0" = "none")) +
  labs(x = "Clone + Maternal Exposure", y = "Infected Individuals",fill="Clone", pattern = "Infection") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
