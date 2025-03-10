library(tidyverse)
library(emmeans)
library(lme4)
library(DHARMa)
source('Bag_data_seq/Bag_ITS_data_prep.R')
library(ggplot2)
library(car)


tax_gs_bag<-read.csv('Processed_data/taxa_w_genome_size_bag_data.csv')

dat_myco_GS<-dat_myco_RA%>%
  left_join(tax_gs_bag%>%select(genus,mean_gs)%>%distinct())%>%
  group_by(Site,Transect,Location)%>%
  summarise(
    CWM_genome = sum(RA_samp * mean_gs, na.rm = TRUE), # Weighted mean genome size
    species_coverage = sum(!is.na(mean_gs)) / n()# Proportion of species with trait data
    ) %>%
  ungroup()%>%
  mutate(log_CWM_genome=log10(CWM_genome+.25/2))%>%
  left_join(dat_myco_RA)

median_coverage <- median(dat_myco_GS$species_coverage, na.rm = TRUE)

dat_myco_GS_filtered <- dat_myco_GS %>%
  filter(species_coverage >= median_coverage)




hist((dat_myco_GS$CWM_genome))
hist((dat_myco_GS$log_CWM_genome))



#fire model
GS_model_fire<-lmer(CWM_genome~ Fire.Severity+Fire.Interval + (1|Site/Transect),  data=dat_myco_GS)

GS_model_fire_filtered<-lmer((CWM_genome)~ Fire.Severity+Fire.Interval + (1|Site/Transect),  data=dat_myco_GS_filtered)

sim_res <- simulateResiduals(fittedModel = GS_model_fire)
plot(sim_res)
#plot(GS_model_fire)
qqPlot(resid(GS_model_fire))
summary(GS_model_fire)
Anova_resin<-round(Anova(GS_model_fire,test='F'), 2) 
Anova_resin
