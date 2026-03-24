
#LOAD LIBRARARIES
library(tidyr)
library(dplyr)
library(vegan)
library(tibble) 
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(fossil)



#From ITS prep script
Soil_Seq_wide<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_Seq_wide.csv')
Soil_data<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Site_Info.csv')
otu_table<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/otu_table_soil.csv',row.names = 'X')

alpha_diversity_ <- otu_table %>%
  as.data.frame() %>%
  rowwise() %>%
  mutate(
    Shannon = diversity(c_across(everything()), index = "shannon"),
    Simpson = diversity(c_across(everything()), index = "simpson"),
    Chao1 = chao1(c_across(everything())),
    Observed = sum(c_across(everything()) > 0),  # Count of non-zero ASVs
    Pielou = Shannon / log(Observed)  # Pielou's evenness
  ) %>%
  ungroup() %>%
  mutate(Pielou = ifelse(Observed <= 1, NA, Pielou))%>%  # Avoid log(1) issue
select(Shannon,Simpson,Chao1,Observed,Pielou)

# Add SampleID back
alpha_diversity<-otu_table%>%
  rownames_to_column(var='barcode')%>%
  mutate(barcode=as.factor(barcode))%>%
  select(barcode) %>% 
  cbind(alpha_diversity_)



alpha_meta_Site_Tran_Loc<-left_join(alpha_diversity,Soil_data)


#Shannon
Shannon<-lmer((Shannon)~  Fire.Severity* Fire.Interval+ (1|Site) , 
                           data=alpha_meta_Site_Tran_Loc)


summary(Shannon)
Anova_Shan<-round(Anova(Shannon,test='F'), 2) 
Anova_Shan
plot(Shannon)
qqPlot(resid(Shannon))

emm_biomass_Interaction<-as.data.frame(emmeans(Simpson,
                                            ~Fire.Severity* Fire.Interval))

emm_biomass_Interaction

Anova_Shan<-Anova_Shan%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Metric='Shannon diversity')

#Pielous
Pielou<-lmer((Pielou)~  Fire.Severity* Fire.Interval+ (1|Site) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Pielou)
Anova_Piel<-round(Anova(Pielou,test='F'), 2) 
Anova_Piel
plot(Pielou)
qqPlot(resid(Pielou))
#r2(Pielou)
emm_biomass_Interval<-as.data.frame(emmeans(Pielou,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Pielou,
                                            ~Fire.Severity))

Anova_Piel<-Anova_Piel%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Metric='Pielou evenness')
#simpson
Simpson<-lmer((Simpson)~  Fire.Severity* Fire.Interval+ (1|Site) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Simpson)
Anova_resin<-round(Anova(Simpson,test='F'), 2) 
Anova_resin
plot(Simpson)
qqPlot(resid(Simpson))
emm_biomass_Interval<-as.data.frame(emmeans(Simpson,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Simpson,
                                            ~Fire.Severity))
#Chao1
Chao1<-lmer(Chao1 ~  Fire.Severity* Fire.Interval+ (1|Site), 
              data=alpha_meta_Site_Tran_Loc)


summary(Chao1)
Anova_Chao<-round(Anova(Chao1,test='F'), 2) 
Anova_Chao
plot(Chao1)
qqPlot(resid(Chao1))
as.data.frame(emmeans(Chao1, ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Chao1,
                                            ~Fire.Severity))
Anova_Chao<-Anova_Chao%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Metric='Chao richness')


Anova_soil<-rbind(Anova_Shan,Anova_Piel,Anova_Chao)%>%
  mutate(source='Soil')
Anova_soil

write.csv(Anova_soil,'Tables/Alpha_diversity_Anova_soil.csv', row.names = FALSE)
