###############################################################################
# Script: Decoupled responses of mycorrhizal fungal communities and function to recurrent wildfire
# Description: Analyses of mycorrhizal biomass, community structure,
#              and trait responses to fire severity and frequency.
# Author: Solomon Maerowitz-McMahan
# Date: [2025-11-12]
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
         starts_with("ITSall")) 

Data_rnd1<- Data_rnd1%>% #seq data
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
         starts_with("ITSall")) 
Data_rnd2<-Data_rnd2%>%
  mutate(
    Biomass_day = if_else(ratio_myco_nonmyco < 0.5, NA_real_, Biomass_day), #remove biomass from rows where seq data suggests nonmyco taxa dominate
    across(c(Biomass_day, Ortho_P_mg_kg, pH), ~log10(.), .names = "{.col}_log"),
    biomass_kg_ha_day = Biomass_day * (1e+03 / 30)
  )

########################################
# 1. Effect of Fire Regime on Nutrients
########################################

# Models: Orthophosphate, Nitrate, Ammonia
m1 <- lmer(Ortho_P_mg_kg_log ~ Fire.Severity + Fire.Interval + (1|Site/Transect), data = Data_rnd1)
m2 <- lmer(Nitrate_mg_kg_log ~ Fire.Severity + Fire.Interval + (1|Site/Transect), data = Data_rnd1)
m3 <- lmer(Ammonia_mg_kg_log ~ Fire.Severity + Fire.Interval + (1|Site/Transect), data = Data_rnd1)

# Combine ANOVA results
Anova_Ortho <- Anova(m1, test = 'F') %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Nutrient = 'Orthophosphate', Round = 'First')
Anova_Nitr  <- Anova(m2, test = 'F') %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Nutrient = 'Nitrate', Round = 'First')
Anova_Ammon <- Anova(m3, test = 'F') %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Nutrient = 'Ammonia', Round = 'First')

Nutrient_1st_Rnd_Anova <- bind_rows(Anova_Ortho, Anova_Nitr, Anova_Ammon)

# Second round resin P
m_Phos <- lmer(log10(Ortho_P_mg_kg) ~ Fire.Interval + (1|Site/Transect), data = Data_rnd2)
Anova_resin_Phos <- Anova(m_Phos, test = 'F') %>%
  round(2) %>% as.data.frame() %>%
  rownames_to_column("Factor") %>%
  mutate(Nutrient = "Orthophosphate", Round = "Second")

# Combined table
Nutrient__Anova <- bind_rows(Nutrient_1st_Rnd_Anova, Anova_resin_Phos)
Nutrient__Anova %>% arrange(Factor)
####################################
# 2. Mycorrhizal Community Analyses
####################################

## 2.1 Alpha Diversity – First Round

# Diversity metrics
alpha_diversity <- Data_rnd1 %>%
  select(Tube_ID, starts_with("ITSall")) %>%
  rowwise() %>%
  mutate(
    Shannon = diversity(c_across(starts_with("ITSall")), index = "shannon"),
    Simpson = diversity(c_across(starts_with("ITSall")), index = "simpson"),
    Chao1 = chao1(c_across(starts_with("ITSall"))),
    Observed = sum(c_across(starts_with("ITSall")) > 0),
    Pielou = if_else(Observed > 1, Shannon / log(Observed), NA_real_)
  ) %>%
  ungroup() %>%
  select(Tube_ID, Shannon, Simpson, Chao1, Observed, Pielou)

# Join metadata
alpha_meta <- alpha_diversity %>%
  left_join(Data_rnd1) %>%
  filter(ratio_myco_nonmyco > 0.5)

# Function to model alpha diversity
analyze_alpha <- function(metric) {
  model <- lmer(as.formula(paste0(metric, " ~ Fire.Severity * Fire.Interval + (1|Site/Transect)")), data = alpha_meta)
  Anova(model, test = "F") %>% round(2) %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Metric = metric)
}

anova_results_1st <- bind_rows(
  analyze_alpha("Shannon"),
  analyze_alpha("Pielou"),
  analyze_alpha("Chao1")
) %>% mutate(source = "1st mycelial collection")

## 2.2 Beta Diversity – First Round

mat_myco <- Data_rnd1 %>% select(starts_with("ITSall"))
permanova_res <- adonis2(mat_myco ~ Fire.Severity * Fire.Interval, data = Data_rnd1, distance = 'robust.aitchison', by = 'margin') %>%
  as.data.frame() %>% rownames_to_column("Factor") %>%
  mutate(Sample_Type = "1st mycelial collection")

## 2.3 CAP Ordination – First Round

cap_mod <- capscale(mat_myco ~ Fire , data = Data_rnd1, distance = 'robust.aitchison', add = TRUE)
anova(cap_mod)
summary(cap_mod)
Cap1_aov<-as.data.frame( anova(cap_mod, by = "margin"))%>%
  rownames_to_column()
cap.allscrs_site <- scores(cap_mod, tidy = TRUE) %>% filter(score == "sites")
prop_var <- round(cap_mod$CCA$eig / cap_mod$tot.chi * 100, 1)





#################################################################################################################################
# Extract scores
scrs <- scores(cap_mod, tidy = TRUE)
scrs_site <- filter(scrs, score == "sites")

# Plot CAP ordination with ellipses
interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_linetypes <- c("High" = "solid", "Low" = "dashed")

cbind(Data_rnd1, scrs_site) %>%
  ggplot(aes(x = CAP1, y = CAP2)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_point(aes(colour = Fire.Interval, shape = Fire.Severity), size = 4, stroke = 1.5) +
  stat_ellipse(aes(colour = Fire.Interval, linetype = Fire.Severity), linewidth = 1) +
  scale_colour_manual(values = interval_colors) +
  scale_linetype_manual(values = severity_linetypes, guide = "none") +
  scale_shape_manual(values = c("High" = 16, "Low" = 1)) +
  labs(
    x = paste0("CAP1 (", prop_var[1], "%)"),
    y = paste0("CAP2 (", prop_var[2], "%)"),
    colour = "Fire Frequency",
    linetype = "Fire Severity",
    shape = "Fire Severity"
  ) +
  theme_classic() +
  theme( legend.position ='top',
    axis.text.x = element_text(hjust = 0.5, size = 12, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    axis.title.x = element_text(size = 14, color = 'black'),
    axis.title.y = element_text(size = 14, color = 'black'),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 10)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 8)),
    colour = guide_legend(override.aes = list(size = 8))
  ) -> p3

p3

ggsave(
  filename = "plots/submission/Figure_2.pdf",
  plot = p3,
  device = "pdf",
  width = 18,          # double-column width
  height = 12,         # proportionate height
  units = "cm"
)

ggsave(
  filename = "plots/submission/Figure_2.png",
  plot = p3,
  device = "png",
  width = 18,          # double-column width
  height = 12,         # proportionate height
  units = "cm"
)


#Mycorrhizal biomasss############

# Compare additive vs interaction model
m_additive <- lmer(Biomass_day_log ~ Fire.Severity + Fire.Interval +
                     Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg  +
                     perc_myco_host_freq + (1 | Site/Transect),
                   data = Data_rnd1)

m_interaction <- lmer(Biomass_day_log ~ Fire.Severity * Fire.Interval +
                        Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg  +
                        perc_myco_host_freq + (1 | Site/Transect),
                      data = Data_rnd1)


# Compare host effect
# m_additive_both_host <- lmer(Biomass_day_log ~ Fire.Severity + Fire.Interval +
#                                Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg  + 
#                                AM_host + EcM_host + (1 | Site/Transect),
#                              data = Data_rnd1) 

# Compare AICs
AIC(m_additive, m_interaction)

# Final selected model
m_biomass_day_resins <- m_additive

# Extract fixed effects
fixed_effects <- tidy(m_biomass_day_resins, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(term=str_remove(term,'Low|Short')) %>% 
  select(term,estimate,std.error)
# Type II ANOVA
anova_df <- Anova(m_biomass_day_resins, test = "F") %>%
  as.data.frame() %>%
  rownames_to_column("Predictor")

anova_df
# Combine slope estimates with ANOVA results
summary_table <- fixed_effects %>%
  rename(Slope = estimate, Std_Error = std.error) %>%
  left_join(anova_df[, c("Predictor", "Df", "Df.res", "F", "Pr(>F)")], by = c("term"="Predictor")) %>%
  rename(p = `Pr(>F)`) %>%
  mutate(
    term = term %>%
      str_replace_all("_mg_kg", "") %>%
      str_replace_all("_", " ") %>%
      str_trim()
  )

# Final cleaned table
summary_table_1st <- summary_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(Round = '1st') %>%
  relocate(Round, term, Slope, Std_Error,F, Df, Df.res, p)
summary_table_1st

Data_rnd1 %>% 
  group_by(Fire.Severity) %>% 
  summarise(
    n = sum(!is.na(biomass_kg_ha_day)),
    all_mean = mean(biomass_kg_ha_day, na.rm = TRUE) * 150,
    se = (sd(biomass_kg_ha_day, na.rm = TRUE) / sqrt(n)) * 150,
    min_biomass = min(biomass_kg_ha_day, na.rm = TRUE)* 150,
    max_biomass = max(biomass_kg_ha_day, na.rm = TRUE)* 150
  )

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
                    Fire.Interval + Fire.Severity + perc_myco_host_freq+ (1 | Site), 
                  data = Bag_data)

summary(C_N_model)
Anova(C_N_model, test = "F")

#----------------#
# C:P Ratio Model
#----------------#
C_P_model <- lmer(log10(C_P) ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                    Fire.Interval + Fire.Severity +  perc_myco_host_freq+ (1 | Site), 
                  data = Bag_data)

summary(C_P_model)
Anova(C_P_model, test = "F")

#----------------#
# N:P Ratio Model
#----------------#
N_P_model <- lmer(log10(N_P) ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                    Fire.Interval + Fire.Severity +  perc_myco_host_freq+(1 | Site), 
                  data = Bag_data)

summary(N_P_model)
Anova(N_P_model, test = "F")

#------------------------------#
# Hyphal Carbon Content Model
#------------------------------#
Carb_Hyph_model <- lmer(Carbon ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                          Fire.Interval + Fire.Severity + perc_myco_host_freq+ (1 | Site), 
                        data = Bag_data)

summary(Carb_Hyph_model)
Anova(Carb_Hyph_model, test = "F")

#-------------------------------#
# Hyphal Nitrogen Content Model
#-------------------------------#
Nitrog_Hyph_model <- lmer(Nitrogen ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                            Fire.Interval + Fire.Severity + perc_myco_host_freq+ (1 | Site), 
                          data = Bag_data)

summary(Nitrog_Hyph_model)
Anova(Nitrog_Hyph_model, test = "F")

#-----------------------------#
# Hyphal Phosphorus Model
#-----------------------------#
Phos_Hyph_model <- lmer(Percent_Phos_ ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                          Fire.Interval + Fire.Severity +  perc_myco_host_freq+ (1 | Site), 
                        data = Bag_data )

summary(Phos_Hyph_model)
Anova(Phos_Hyph_model, test = "F")

####2nd mycelial Collection##############

# Fit LMM
m_biomass <- lmer(Biomass_day_log ~ Fire.Interval + Ortho_P_mg_kg + perc_myco_host_freq+
                    (1 | Site / Transect), data = Data_rnd2)

# Summary diagnostics
summary(m_biomass)
Anova_biomass <- Anova(m_biomass, test = "F") %>% round(2)
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

# Save summary table
#write_csv(summary_table, "Tables/Biomass_Anova_2nd_rnd.csv")


###Community###
#Second mycelial collection###

#alpha diversity######

# Extract community matrix
otu_table <- Data_rnd2 %>% select(starts_with("ITSall"))

# Calculate alpha diversity
alpha_div <- Data_rnd2 %>%
  select(ID, starts_with("ITSall")) %>%
  rowwise() %>%
  mutate(
    Shannon = diversity(c_across(starts_with("ITSall")), index = "shannon"),
    Simpson = diversity(c_across(starts_with("ITSall")), index = "simpson"),
    Chao1 = chao1(c_across(starts_with("ITSall"))),
    Observed = sum(c_across(starts_with("ITSall")) > 0),
    Pielou = if_else(Observed > 1, Shannon / log(Observed), NA_real_)
  ) %>%
  ungroup() %>%
  select(ID, Shannon, Simpson, Chao1, Observed, Pielou)

# Join metadata
alpha_meta <- left_join(alpha_div, Data_rnd2) %>%
  select(ID, Shannon, Simpson, Chao1, Observed, Pielou,
         Site, Transect, Location, Fire.Interval)

# Function to model and return ANOVA table
run_alpha_model <- function(metric, data) {
  formula <- as.formula(paste0(metric, " ~ Fire.Interval + (1|Site/Transect)"))
  model <- lmer(formula, data = data)
  print(summary(model))
  print(qqPlot(resid(model), main = paste("QQ plot:", metric)))
  anova_res <- Anova(model, test = "F") %>%
    round(2) %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(Metric = metric)
  return(anova_res)
}

# Run models and combine results
anova_results <- bind_rows(
  run_alpha_model("Shannon", alpha_meta),
  run_alpha_model("Pielou", alpha_meta),
  run_alpha_model("Chao1", alpha_meta)
) %>%
  mutate(source = "2nd mycelial collection")

anova_results

#Beta diversity#########

# Filter out samples with no mycorrhizal reads
Data_rnd2 <- Data_rnd2 %>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

# Extract community matrix
mat_myco <- Data_rnd2 %>% select(starts_with("ITSall"))

# PERMANOVA
adonis_rnd2 <- adonis2(mat_myco ~ Fire.Interval, data = Data_rnd2,
                      distance = "robust.aitchison", by = "margin") %>%
  as.data.frame() %>%
  rownames_to_column("Factor") %>%
  mutate(Sample_Type = "2nd mycelial collection")

adonis_rnd2
# CAP Ordination
cap_model <- capscale(mat_myco ~ Fire.Interval, data = Data_rnd2,
                      distance = "robust.aitchison", add = TRUE)

anova(cap_model)
summary(cap_model)

# Axis proportions
proportions <- round(cap_model$CCA$eig / cap_model$tot.chi * 100, 1)

# Extract tidy scores
scrs <- scores(cap_model, tidy = TRUE)

scrs_site <- scrs %>% filter(score == "sites")
scrs_cent <- scrs %>% filter(score == "centroids")

# Define colors
interval_colors <- c("Long" = "darkred", "Short" = "orange")

# Plot CAP
p2 <- cbind(Data_rnd2, scrs_site) %>%
  ggplot(aes(x = CAP1, y = MDS1 )) +
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_point(aes(colour = Fire.Interval), size = 8) +
  stat_ellipse(aes(color = Fire.Interval), level = 0.95, linewidth = 2) +
  scale_colour_manual(values = interval_colors) +
  labs(
    x = paste0("CAP1 (", proportions[1], "%)"),
    y = "MDS1",
    colour = "Fire frequency"
  ) +
  xlim(min(scrs_site$CAP1) - 1, max(scrs_site$CAP1) + 0.5) +
  theme_bw(base_size = 20) +
  theme(
    axis.text = element_text(color = "black", size = 36),
    axis.title = element_text(color = "black", size = 36),
    legend.position = "top",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 8)))

p2







