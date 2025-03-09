library(tidyverse)
library(emmeans)
library(lme4)
library(DHARMa)
source('Bag_data_seq/Bag_ITS_data_prep.R')
library(ggplot2)
library(car)
library(glmmTMB)


dat_myco_RA<-dat_myco_RA%>%
  mutate(log_RA_total=log10(RA_total_reads))
#this is relative abundance compared to other mycorrhizal fungi



hist(log10(dat_myco_RA$RA_total_reads))



#fire model
RA_model<-lmer(log_RA_total~ Fire.Severity+Fire.Interval + 
                         Nitrate_mg_kg+Ammonia_mg_kg+Ortho_P_mg_kg +
                         perc_myco_host_freq+
                         (1|Site/Transect/Location),
                       data=dat_myco_RA)

RA_model <- glmmTMB(
  RA_total_reads ~ Fire.Severity + Fire.Interval + 
    Nitrate_mg_kg + Ammonia_mg_kg + Ortho_P_mg_kg +
    perc_myco_host_freq +
    (1 | Site/Transect/Location), 
  data = dat_myco_RA, 
  family = tweedie)  # Models abundance given presence


sim_res <- simulateResiduals(fittedModel = RA_model)
plot(sim_res)
#plot(RA_model)
qqPlot(resid(RA_model))
summary(RA_model)
Anova_resin<-round(Anova(RA_model,test='F'), 2) 
Anova_resin