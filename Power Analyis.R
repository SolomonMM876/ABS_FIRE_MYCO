

library(pwr)
library(simr)
library(effectsize)

Bag_Site

m1<-lmer(myc_bag_yield_est~Severity+ Interval + (1|Site) , data=Bag_Site)
summary(m1)




# Calculate effect sizes for fixed effects
effect_sizes <- effectsize::standardize_parameters(m1)
print(effect_sizes)

m1_sim <- extend(m1, along="Site", n=30) # Assuming extending to 30 sites for power simulation

# Power analysis for fixed effects
power_severity <- powerSim(m1_sim, fixed("Severity"), nsim=1000)
power_interval <- powerSim(m1_sim, fixed("Interval"), nsim=1000)

# Print the power analysis results
print(power_severity)
print(power_interval)

var_comp <- as.data.frame(VarCorr(m1))
print(var_comp)




# Extend the model by increasing the number of samples per site
# Here, we double the number of samples per site from 8 to 16
new_n <- 64 / 8  # Doubling the number of samples

m1_extended <- extend(m1, within = "Site", n = new_n * length(unique(Bag_Site$Interval)) * length(unique(Bag_Site$Severity)))

# Note: Adjust the `new_n` factor based on how much you want to increase the samples per site.


# Perform power analysis for the Interval predictor in the extended model
power_result_interval_extended <- powerSim(m1_extended, fixed("Interval"), nsim = 500)

# Perform power analysis for the Severity predictor in the extended model
power_result_severity_extended <- powerSim(m1_extended, fixed("Severity"), nsim = 500)

# Print the power analysis result for the Interval predictor
print(power_result_interval_extended)

# Print the power analysis result for the Severity predictor
print(power_result_severity_extended)










# Fit the ANOVA model with the specified data and factors
anova_model <- aov(myc_bag_yield_est  ~ Interval + Severity, data = Bag_Site)
anova_summary <- summary(anova_model)
anova_table <- anova_summary[[1]]

# Calculate partial eta squared for each effect
SS_effects <- anova_table$`Sum Sq`[-length(anova_table$`Sum Sq`)] # Exclude residuals
SS_error <- anova_table$`Sum Sq`[length(anova_table$`Sum Sq`)]
eta_squared_partial <- SS_effects / (SS_effects + SS_error)

# Print partial eta squared values for each effect
print(eta_squared_partial)

# Convert partial eta squared to Cohen's f for each effect
f_effects <- sqrt(eta_squared_partial / (1 - eta_squared_partial))

# Number of levels for each factor
levels_interval <- length(unique(Bag_Site$Interval))
levels_severity <- length(unique(Bag_Site$Severity))

# Calculate sample size per cell
n_per_cell <- length(Bag_Site$myc_bag_yield_est) / (levels_interval * levels_severity)

# Perform post hoc power analysis for each effect
# For main effect of Interval
power_result_interval <- pwr.anova.test(k = levels_interval, 
                                        n = n_per_cell, 
                                        f = f_effects[1], 
                                        sig.level = 0.05)

# For main effect of Severity
power_result_severity <- pwr.anova.test(k = levels_severity, 
                                        n = n_per_cell, 
                                        f = f_effects[2], 
                                        sig.level = 0.05)

# For interaction between Interval and Severity
power_result_interaction <- pwr.anova.test(k = levels_interval * levels_severity, 
                                           n = n_per_cell, 
                                           f = f_effects[3], 
                                           sig.level = 0.05)

# Print the power analysis results
print(power_result_interval)
print(power_result_severity)
print(power_result_interaction)

