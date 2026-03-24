library(lme4)
library(simr)
library(effectsize)
library(pwr)

# Load Data
Bag_data <- read.csv('Processed_data/All_Bag_data.csv')

# Fit Linear Mixed Model
m_biomass <- lmer(log10_biomass_day ~ Fire.Severity + Fire.Interval +
                    Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + pH + 
                    perc_myco_host_freq +
                    (1|Site/Transect), 
                  data = Bag_data)

# Calculate effect sizes for fixed effects
effect_sizes <- effectsize::standardize_parameters(m_biomass)
print(effect_sizes)

# Check variance components
var_comp <- as.data.frame(VarCorr(m_biomass))
print(var_comp)

# Baseline Power Analysis for Fire Severity and Fire Interval
power_severity <- powerSim(m_biomass, fixed("Fire.Severity"), nsim = 1000)
power_interval <- powerSim(m_biomass, fixed("Fire.Interval"), nsim = 1000)

print(power_severity)
print(power_interval)

# Simulate adding more sites to estimate required sample size
# Extend to 30, 40, and 50 sites for comparison
m1_30 <- extend(m_biomass, along = "Site", n = 30)
m1_40 <- extend(m_biomass, along = "Site", n = 40)
m1_60 <- extend(m_biomass, along = "Site", n = 60)

# Power Analysis for extended models
power_30 <- powerSim(m1_30, fixed("Fire.Interval"), nsim = 500)
power_40 <- powerSim(m1_40, fixed("Fire.Interval"), nsim = 500)


power_60 <- powerSim(m1_60, fixed("Fire.Severity"), nsim = 500)

print(power_30)
print(power_40)
print(power_50)

# Extend within-site sampling: double samples per site
n_per_site <- nrow(Bag_data) / length(unique(Bag_data$Site))
m1_extended <- extend(m_biomass, within = "Site", n = n_per_site * 2)

power_extended <- powerSim(m1_extended, fixed("Fire.Severity"), nsim = 500)
print(power_extended)

# -------------------------------------------------
# ANOVA-Based Power Analysis
# -------------------------------------------------

# Fit ANOVA model
anova_model <- aov(myc_bag_yield_est ~ Fire.Interval + Fire.Severity, data = Bag_data)
anova_summary <- summary(anova_model)
anova_table <- anova_summary[[1]]

# Partial eta squared
SS_effects <- anova_table$`Sum Sq`[-length(anova_table$`Sum Sq`)] # Exclude residuals
SS_error <- anova_table$`Sum Sq`[length(anova_table$`Sum Sq`)]
eta_squared_partial <- SS_effects / (SS_effects + SS_error)

print(eta_squared_partial)

# Convert to Cohen’s f
f_effects <- sqrt(eta_squared_partial / (1 - eta_squared_partial))

# Number of levels per factor
levels_interval <- length(unique(Bag_data$Fire.Interval))
levels_severity <- length(unique(Bag_data$Fire.Severity))

# Sample size per cell
n_per_cell <- length(Bag_data$myc_bag_yield_est) / (levels_interval * levels_severity)

# Power analysis for each effect
power_result_interval <- pwr.anova.test(k = levels_interval, 
                                        n = n_per_cell, 
                                        f = f_effects[1], 
                                        sig.level = 0.05)

power_result_severity <- pwr.anova.test(k = levels_severity, 
                                        n = n_per_cell, 
                                        f = f_effects[2], 
                                        sig.level = 0.05)

power_result_interaction <- pwr.anova.test(k = levels_interval * levels_severity, 
                                           n = n_per_cell, 
                                           f = f_effects[3], 
                                           sig.level = 0.05)

print(power_result_interval)
print(power_result_severity)
print(power_result_interaction)

