library(tidyverse)
library(emmeans)
library(lme4)
library(DHARMa)
source('Bag_data_seq/Bag_ITS_data_prep.R')
library(ggplot2)
library(car)


dat_explo<-dat_myco_RA%>%
  group_by(Site,Transect,Location,exploration_type)%>%
  summarise(explo_count=sum(count))%>%
  left_join(dat_myco_RA%>%select(Site,Transect,Location,Fire.Severity,Fire.Interval,reads_samp)%>%distinct())%>%
  left_join(Bag_data%>%select(Site,Transect,Location,pH:Ortho_P_mg_kg))%>%
  mutate(RA_explo=explo_count/reads_samp,
         log_RA_explo=log10(RA_explo))
#this is relative abundance compared to other mycorrhizal fungi



hist(log10(dat_explo$RA_explo))



#fire model
explo_model_fire<-lmer((log_RA_explo)~ exploration_type *(Fire.Severity+Fire.Interval) + 
                         (1|Site/Transect),
                       data=dat_explo)



sim_res <- simulateResiduals(fittedModel = explo_model_fire)
plot(sim_res)
#plot(explo_model_fire)
qqPlot(resid(explo_model_fire))
summary(explo_model_fire)
Anova_resin<-round(Anova(explo_model_fire,test='F'), 2) 
Anova_resin

#nutrient model
explo_model_nutri<-lmer((log_RA_explo)~ exploration_type *(Nitrate_mg_kg+Ammonia_mg_kg+Ortho_P_mg_kg) + 
                         (1|Site/Transect),
                       data=dat_explo)



sim_res <- simulateResiduals(fittedModel = explo_model_nutri)
plot(sim_res)
qqPlot(resid(explo_model_nutri))
summary(explo_model_nutri)
Anova_resin<-round(Anova(explo_model_nutri,test='F'), 2) 
Anova_resin

