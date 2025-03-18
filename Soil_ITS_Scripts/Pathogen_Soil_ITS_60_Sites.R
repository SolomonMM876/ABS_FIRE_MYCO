source('Soil_ITS_Scripts/Soil_ITS_data prep.R')
library(ggplot2)





dat_path<-dat_all_RA%>%
  rownames_to_column("SH_ID_OTU")%>%
  #keep only OTUs with pathogen
  filter(guild2=='Pathogen')%>%
  group_by(Site,Transect)%>%
  summarise(path_count=sum(readcount))%>%
  left_join(dat_all_RA%>%select(Site,Transect,Severity,Interval,reads_samp)%>%distinct())%>%
  mutate(RA_path=path_count/reads_samp,
        log_RA_path=log10(RA_path+0.0002853067/2))

dat_plant_path<-dat_all_RA%>%
  rownames_to_column("SH_ID_OTU")%>%
  #keep only OTUs with pathogen
  filter(str_detect(guild, 'Plant Pathogen'))%>%
  group_by(Site,Transect)%>%
  summarise(plant_path_count=sum(readcount))%>%
  left_join(dat_all_RA%>%select(Site,Transect,Severity,Interval,reads_samp)%>%distinct())%>%
  mutate(RA_plant_path=plant_path_count/reads_samp,
         log_RA_plant_path=log10(RA_plant_path+0.002353585/2))

dat_arma<-dat_all_RA%>%
  rownames_to_column("SH_ID_OTU")%>%
  #keep only OTUs with pathogen
  filter(Genus=='Armillaria')


hist(dat_path$log_RA_path)


#library(glmmTMB)
library(DHARMa)

#fire model for all paths
path_model_fire<-lmer((log_RA_path)~ (Severity+Interval) + 
                         (1|Site),
                       data=dat_path)




sim_res <- simulateResiduals(fittedModel = path_model_fire)
plot(sim_res)
#plot(explo_model_fire)
qqPlot(resid(path_model_fire))
summary(path_model_fire)
Anova_resin<-round(Anova(path_model_fire,test='F'), 2) 
Anova_resin


emm_severity <- as.data.frame(emmeans(path_model_fire, ~ Severity))
emm_interval <-as.data.frame(emmeans(path_model_fire, ~ Interval))

path_severity<-emm_severity%>%
  ggplot( aes(x =Severity, y = (10^(emmean))) )+
  geom_col(size=3, alpha=.6)+
  geom_errorbar(aes(ymin=(10^(emmean-SE)),ymax=(10^(emmean+SE))), width=.2) +
  labs(x = "Fire Severity ", y = "RA_pathogen") +
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 1, size = 15, face = "bold",  angle = 45),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

path_severity

path_interval<-emm_interval%>%
  ggplot( aes(x =Interval, y = (10^(emmean))) )+
  geom_col(size=3, alpha=.6)+
  geom_errorbar(aes(ymin=(10^(emmean-SE)),ymax=(10^(emmean+SE))), width=.2) +
  labs(x = "Fire Frequency ", y = "RA_pathogen") +
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 1, size = 15, face = "bold",  angle = 45),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

path_interval



#fire model for plant paths

plant_path_model_fire<-lmer((log_RA_plant_path)~ (Severity+Interval) + 
                        (1|Site),
                      data=dat_plant_path)



sim_res <- simulateResiduals(fittedModel = path_model_fire)
plot(sim_res)
#plot(explo_model_fire)
qqPlot(resid(path_model_fire))
summary(path_model_fire)
Anova_resin<-round(Anova(path_model_fire,test='F'), 2) 
Anova_resin


emm_severity <- as.data.frame(emmeans(path_model_fire, ~ Severity))
emm_interval <-as.data.frame(emmeans(path_model_fire, ~ Interval))

path_severity<-emm_severity%>%
  ggplot( aes(x =Severity, y = (10^(emmean))) )+
  geom_col(size=3, alpha=.6)+
  geom_errorbar(aes(ymin=(10^(emmean-SE)),ymax=(10^(emmean+SE))), width=.2) +
  labs(x = "Fire Severity ", y = "RA_pathogen") +
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 1, size = 15, face = "bold",  angle = 45),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

path_severity

path_interval<-emm_interval%>%
  ggplot( aes(x =Interval, y = (10^(emmean))) )+
  geom_col(size=3, alpha=.6)+
  geom_errorbar(aes(ymin=(10^(emmean-SE)),ymax=(10^(emmean+SE))), width=.2) +
  labs(x = "Fire Frequency ", y = "RA_pathogen") +
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 1, size = 15, face = "bold",  angle = 45),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

path_interval
