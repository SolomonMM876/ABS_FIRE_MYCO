library(lme4)
library(ggplot2)
library(car)
library(emmeans)
library(tidyverse)


#First Round
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')



#log transform response variables
Bag_data<-Bag_data%>%
  mutate(Ortho_P_mg_kg_log=log10(Ortho_P_mg_kg),
         Nitrate_mg_kg_log=log10(Nitrate_mg_kg),
         Ammonia_mg_kg_log=log10(Ammonia_mg_kg),
         pH_log=log10(pH))





#Orthophosphate
m_additive<-lmer(log10(Ortho_P_mg_kg)~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_data)

m_interaction<-lmer(log10(Ortho_P_mg_kg)~Fire.Severity* Fire.Interval + (1|Site/Transect) , data=Bag_data)

AIC(m_additive, m_interaction)

m1<-lmer(log10(Ortho_P_mg_kg)~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_data)

#1
summary(m1)
Anova_Ortho<-Anova(m1,test='F')
Anova_Ortho<-Anova_Ortho %>% 
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Nutrient='Orthophosphate',
         Round= 'First')
plot(m1)
qqPlot(resid(m1))
pairs(emmeans(m1, ~ Fire.Interval))
emm_Interval<-as.data.frame(emmeans(m1, ~Fire.Interval))
emmip(m1, ~Fire.Interval)

#g/hectare 10,000 m^2 ×0.1 m= 1,000 m^3 
# 1,000m^3 x (x mg /15cm^3) x (1g/1000mg) x 1000kg/ton= kg/ha
(1000/15)


#Nitrate
m_additive<-lmer(log10(Nitrate_mg_kg)~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_data)

m_interaction<-lmer(log10(Nitrate_mg_kg)~Fire.Severity* Fire.Interval + (1|Site/Transect) , data=Bag_data)

AIC(m_additive, m_interaction)

m2<-lmer(log10(Nitrate_mg_kg)~  Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_data)
#perc_N_fix_freq

#m2<-lmer((Nitrate_mg_kg_log)~Fire.Severity+ Fire.Interval + Ortho_P_mg_kg + Ammonia_mg_kg + pH+ (1|Site/Transect) , data=Bag_data)


#1
summary(m2)
Anova_Nitr<- Anova(m2,test='F')
Anova_Nitr<-Anova_Nitr %>% 
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Nutrient='Nitrate',
         Round= 'First')
plot(m2)
qqPlot(resid(m2))
pairs(emmeans(m2, ~ Fire.Interval))
emm_Interval<-as.data.frame(emmeans(m2, ~Fire.Interval))
emm_Interval <- emm_Interval %>%
  mutate(
    emmean_bt = 10^emmean,
    lower.CL_bt = 10^lower.CL,
    upper.CL_bt = 10^upper.CL
  )
emmip(m2, ~Fire.Interval)


#g/hectare:
#10,000 m^2 ×0.1 m= 1,000 m^3 
#1,000m^3 x (x mg /15cm^3) x (1g/1000mg) x 1000kg/ton x 1000g/kg= g/ha
(1e+06/15)





#Ammonia
m_additive<-lmer(log10(Ammonia_mg_kg)~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_data)

m_interaction<-lmer(log10(Ammonia_mg_kg)~Fire.Severity* Fire.Interval + (1|Site/Transect) , data=Bag_data)

AIC(m_additive, m_interaction)

m3<-lmer(log10(Ammonia_mg_kg) ~   Fire.Interval * Fire.Severity  + (1|Site/Transect) , data=Bag_data)
#perc_N_fix_freq

#m3<-lmer((Ammonia_mg_kg_log)~Fire.Severity+ Fire.Interval + Ortho_P_mg_kg + Nitrate_mg_kg + pH+ (1|Site/Transect) , data=Bag_data)

#3 (NH4~Regime) Not Sig
summary(m3)
Anova_Ammon<-Anova(m3,test='F')

Anova_Ammon<-Anova_Ammon %>% 
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Nutrient='Ammonia',
         Round= 'First')
plot(m3)
qqPlot(resid(m3))
emm_Interval<-as.data.frame(emmeans(m3, ~Fire.Interval))
emm_Interval <- emm_Interval %>%
  mutate(
    emmean_bt = 10^emmean,
    lower.CL_bt = 10^lower.CL,
    upper.CL_bt = 10^upper.CL
  )


Nutrient_1st_Rnd_Anova<-rbind(Anova_Ortho,Anova_Nitr,Anova_Ammon)

saveRDS(Nutrient_1st_Rnd_Anova,'Tables/Nutients_1st_Rnd_Anova.rds')


############Second Round###############################

combined_data<-read.csv('ABS_Second_Rnd/processed_data/combined_data.csv')


#Resins related to biomass

m_Ammonia<-lmer(log10(Amm_mg_kg) ~  Fire.Interval + Fire.Severity  + (1|Site/Transect) , data=combined_data)


summary(m_Ammonia)
Anova_resin_Amm<-round(Anova(m_Ammonia,test='F'), 2) 
Anova_resin_Amm<-Anova_resin_Amm%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Nutrient= "Ammonia",
    Round='Second')
plot(m_Ammonia)
qqPlot(resid(m_Ammonia))


m_Phos<-lmer(log10(Ortho_P_mg_kg) ~   Fire.Interval + Fire.Severity  + (1|Site/Transect) , data=combined_data)



summary(m_Phos)
Anova_resin_Phos<-round(Anova(m_Phos,test='F'), 2) 
Anova_resin_Phos<-Anova_resin_Phos%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Nutrient= "Orthophosphate",
         Round='Second')
plot(m_Phos)
qqPlot(resid(m_Phos))

Nutrient_2nd_Rnd_Anova<-rbind(Anova_resin_Phos,Anova_resin_Amm)

saveRDS(Nutrient_2nd_Rnd_Anova,'Tables/Nutients_2nd_Rnd_Anova.rds')
