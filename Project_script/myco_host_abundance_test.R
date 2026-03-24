###############################################################################
# Script: Mycorrhizal Community Response to Fire Regime
# Description: Analyses of mycorrhizal biomass, community structure,
#              and trait responses to fire severity and frequency.
# Author: [Your Name]
# Date: [YYYY-MM-DD]
###############################################################################

########################
# Load Required Packages
########################

library(tidyverse)
library(lme4)
library(car)
library(vegan)
library(emmeans)
library(fossil)
library(broom.mixed)

########################
# Load and Prepare Data
########################
#load meta data together
source('df_joins.R')
#join meta data with 1 rnd seq data
source('Bag_data_seq/Bag_ITS_data_prep.R')
rm(list = ls())

# First round data (1st mycelial collection)
Bag_Seq_wide <- read.csv('Processed_data/Bag_Seq_wide.csv')
myco_tax <- read.csv('Processed_data/Bag_dat_myco_tax.csv')

Data_rnd1 <- Bag_Seq_wide %>%
  select(Site, Transect, Location, Fire.Interval, Fire.Severity, #Site
         Latitude, Longitude,#location
         Ortho_P_mg_kg, Nitrate_mg_kg, Ammonia_mg_kg, pH, #soil
         perc_myco_host_freq,AM_host=VAM,EcM_host=ECM,#myco veg
         Biomass_day, Tube_ID, #biomass
         ratio_myco_nonmyco, myco_reads,#myco data
         starts_with("ITSall")) %>% #seq data
  mutate(
    Biomass_day = if_else(ratio_myco_nonmyco < 0.5, NA_real_, Biomass_day), #remove biomass from rows where seq data suggests nonmyco taxa dominate
    across(c(Biomass_day, Ortho_P_mg_kg, Nitrate_mg_kg, Ammonia_mg_kg, pH), ~log10(.), .names = "{.col}_log"),
    # Convert to kg/ha/day (10,000 m^2 × 0.1 m) x (mg hyph /(15cm^3)x 2 bags ) x (1g/1000mg) x 1000kg/ton = kg/ha
    biomass_kg_ha_day = Biomass_day * (1e+03 / 30),
    Fire = paste(Fire.Severity, "x", Fire.Interval, sep = "\n")
  )
# Second round data (2nd mycelial collection)
Bag_Seq_wide <- read_csv('ABS_Second_Rnd/Processed_data/Bag_Seq_wide.csv')
myco_tax <- read_csv('ABS_Second_Rnd/Processed_data/Bag_dat_myco_tax.csv')

Data_rnd2 <- Bag_Seq_wide %>%
  select(Site, Transect, Location, ID, Fire.Interval,
         Ortho_P_mg_kg, pH = Avg_pH,
         perc_myco_host_freq,AM_host=VAM,EcM_host=ECM,#myco veg
         Biomass_day,
         ratio_myco_nonmyco, myco_reads,
         starts_with("ITSall")) %>%
  mutate(
    Biomass_day = if_else(ratio_myco_nonmyco < 0.5, NA_real_, Biomass_day), #remove biomass from rows where seq data suggests nonmyco taxa dominate
    across(c(Biomass_day, Ortho_P_mg_kg, pH), ~log10(.), .names = "{.col}_log"),
    biomass_kg_ha_day = Biomass_day * (1e+03 / 30)
  )


#Mycorrhizal biomasss############


# Compare additive vs interaction model
# m_additive <- lmer(Biomass_day_log ~ Fire.Severity + Fire.Interval +
#                      Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg  +
#                      perc_myco_host_freq + (1 | Site/Transect),
#                    data = Data_rnd1)


# Compare host effect
m_additive_both_host <- lmer(Biomass_day_log ~ Fire.Severity + Fire.Interval +
                              Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg  + 
                              AM_host + EcM_host + (1 | Site/Transect),
                            data = Data_rnd1)

# Extract fixed effects
fixed_effects <- tidy(m_additive_both_host, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(term=str_remove(term,'Low|Short')) %>% 
  select(term,estimate,std.error)
# Type II ANOVA
 Anova(m_additive_both_host, test = "F") 

# Fit LMM
m_biomass <- lmer(Biomass_day_log ~ Fire.Interval + Ortho_P_mg_kg +
                    AM_host + EcM_host +
                    (1 | Site / Transect), data = Data_rnd2)

# Summary diagnostics
summary(m_biomass)
 Anova(m_biomass, test = "F") %>% round(2)
qqPlot(resid(m_biomass), main = "QQ Plot: Biomass Model")
plot(m_biomass)

# Tidy fixed effects and ANOVA
fixed_effects <- tidy(m_biomass, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(Predictor = str_remove_all(term, "Low|Short"),
         Predictor = str_replace_all(Predictor, "_mg_kg", ""),
         Predictor = str_replace_all(Predictor, "_", " ")) %>%
  select(Predictor, estimate, std.error)

anova_df <- Anova_biomass %>%
  as.data.frame() %>%
  rownames_to_column("Predictor") %>%
  mutate(Predictor = str_replace_all(Predictor, "_mg_kg", ""),
         Predictor = str_replace_all(Predictor, "_", " "),
         Predictor = trimws(Predictor))

# Combine slope estimates with ANOVA results
summary_table <- fixed_effects %>%
  rename(Slope = estimate, Std_Error = std.error) %>%
  left_join(anova_df, by = "Predictor") %>%
  rename(p = `Pr(>F)`) %>%
  mutate(across(where(is.numeric), round, 2),
         Round = "2nd") %>%
  relocate(Round, Predictor, F, Slope, Std_Error, Df, Df.res, p)




#Traits of mycorrhizal fungi ###########

# Load and prepare data
CNP_clean <- read.csv("Processed_data/CNP_clean.csv")

Bag_data <- Data_rnd1 %>%
  group_by(Site, Transect) %>%
  summarise(across(c(Ortho_P_mg_kg, Nitrate_mg_kg, Ammonia_mg_kg), ~mean(.x, na.rm = TRUE))) %>%
  left_join(CNP_clean, by = c("Site", "Transect")) %>% 
  left_join(Data_rnd1 %>% distinct(Site,Transect,AM_host,EcM_host,perc_myco_host_freq))

#----------------#
# C:N Ratio Model
#----------------#
C_N_model <- lmer(log10(C_N) ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                    Fire.Interval + Fire.Severity +AM_host + EcM_host+ (1 | Site), 
                  data = Bag_data)

summary(C_N_model)
Anova(C_N_model, test = "F")

#----------------#
# C:P Ratio Model
#----------------#
C_P_model <- lmer(log10(C_P) ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                    Fire.Interval + Fire.Severity +AM_host + EcM_host+ (1 | Site), 
                  data = Bag_data)

summary(C_P_model)
Anova(C_P_model, test = "F")

#----------------#
# N:P Ratio Model
#----------------#
N_P_model <- lmer(log10(N_P) ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                    Fire.Interval + Fire.Severity +AM_host + EcM_host+ (1 | Site), 
                  data = Bag_data)

summary(N_P_model)
Anova(N_P_model, test = "F")

#------------------------------#
# Hyphal Carbon Content Model
#------------------------------#
hist(Bag_data$Carbon)
Carb_Hyph_model <- lmer(Carbon ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                          Fire.Interval + Fire.Severity +AM_host + EcM_host+ (1 | Site), 
                        data = Bag_data)

summary(Carb_Hyph_model)
Anova(Carb_Hyph_model, test = "F")

#-------------------------------#
# Hyphal Nitrogen Content Model
#-------------------------------#
Nitrog_Hyph_model <- lmer(Nitrogen ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                            Fire.Interval + Fire.Severity +AM_host + EcM_host+ (1 | Site), 
                          data = Bag_data)

summary(Nitrog_Hyph_model)
Anova(Nitrog_Hyph_model, test = "F")

#-----------------------------#
# Hyphal Phosphorus Model
#-----------------------------#
Phos_Hyph_model <- lmer(Percent_Phos_ ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                          Fire.Interval + Fire.Severity +AM_host + EcM_host+ (1 | Site), 
                        data = Bag_data )

summary(Phos_Hyph_model)
Anova(Phos_Hyph_model, test = "F")



