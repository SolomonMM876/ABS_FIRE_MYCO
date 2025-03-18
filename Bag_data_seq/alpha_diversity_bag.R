
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

#source('Bag_data_seq/Bag_ITS_data_prep.R')


#All meta data from 12 sites with bags collected
#From ITS prep script
Bag_Seq_wide<-read.csv('Processed_data/Bag_Seq_wide.csv')
Bag_data<-read.csv('Processed_data/Updated_Bag_data.csv')
dat_summary <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_summarised.tsv')
otu_table<-read.csv('Processed_data/otu_table_bag.csv',row.names = 'X')

alpha_diversity <- otu_table %>%
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
alpha_diversity<-alpha_diversity%>%
  rownames_to_column(var='Tube_ID')%>%
  mutate(Tube_ID=as.numeric(Tube_ID))


# alpha_diversity<-alpha_diversity%>%
#   filter(!(Shannon  == 0 & Simpson  == 0))

# tmp<-Bag_data%>%
#   select(Tube_ID,Fire.Interval,Fire.Severity,Site, Transect, Location)

#anti_join(tmp,wide_seq)

alpha_meta_Site_Tran_Loc<-left_join(alpha_diversity,Bag_data)

#alpha_meta_Site<-left_join(Bag_data%>%select(-Transect,-Location), alpha_diversity)

# alpha_meta_Site_Tran_Loc<-alpha_meta_Site_Tran_Loc%>%
#   filter( !is.na(Shannon|Simpson))

#Shannon
Shannon<-lmer((Shannon)~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                Nitrate_mg_kg+ Ammonia_mg_kg+ perc_myco_host_freq+ perc_N_fix_freq+ (1|Site/Transect) , 
                           data=alpha_meta_Site_Tran_Loc)


summary(Shannon)
Anova_resin<-round(Anova(Shannon,test='F'), 2) 
Anova_resin
plot(Shannon)
qqPlot(resid(Shannon))
r2(Shannon)

#Pielous
Pielou<-lmer((Pielou)~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                Nitrate_mg_kg+ Ammonia_mg_kg+perc_myco_host_freq+ perc_N_fix_freq+(1|Site/Transect) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Pielou)
Anova_resin<-round(Anova(Pielou,test='F'), 2) 
Anova_resin
plot(Pielou)
qqPlot(resid(Pielou))
#r2(Pielou)
emm_biomass_Interval<-as.data.frame(emmeans(Pielou,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Pielou,
                                            ~Fire.Severity))
#simpson
Simpson<-lmer((Simpson)~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                Nitrate_mg_kg+ Ammonia_mg_kg+ (1|Site/Transect) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Simpson)
Anova_resin<-round(Anova(Simpson,test='F'), 2) 
Anova_resin
plot(Simpson)
qqPlot(resid(Simpson))
r2(Simpson)
emm_biomass_Interval<-as.data.frame(emmeans(Simpson,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Simpson,
                                            ~Fire.Severity))
#Chao1
Chao1<-lmer(Chao1 ~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                 Nitrate_mg_kg+ Ammonia_mg_kg+ perc_myco_host_freq+ (1|Site/Transect) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Chao1)
Anova_resin<-round(Anova(Chao1,test='F'), 2) 
Anova_resin
plot(Chao1)
qqPlot(resid(Chao1))
r2(Chao1)
emm_biomass_Interval<-as.data.frame(emmeans(Chao1,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Chao1,
                                            ~Fire.Severity))
####Plots########

library(ggplot2)

ggplot(alpha_meta_Site_Tran_Loc, aes(x = perc_myco_host_freq, y = Shannon)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatterplot points
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Linear regression line
  labs(x = "Mycorrhizal Host Frequency (%)", y = "Shannon Diversity", 
       title = "Shannon Diversity vs. Mycorrhizal Host Frequency") +
  theme_minimal()
