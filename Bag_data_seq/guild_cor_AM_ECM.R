source('Bag_data_seq/Bag_ITS_data_prep.R')


# Summarize abundance of each guild per sample
guild_abundance <- dat_myco_RA_bag %>%
  group_by(Site,Transect,Location, Fire.Interval, Fire.Severity, guild2) %>%
  summarise(total_abundance = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = guild2, values_from = total_abundance, values_fill = 0)%>%
  rename(AM=`Arbuscular Mycorrhizal`)

# Rename columns for clarity

library(rstatix)

# Compute correlations per Fire Interval
cor_results <- guild_abundance %>%
  group_by(Fire.Interval) %>%
  summarise(correlation = cor(Ectomycorrhizal, AM, method = "pearson", use = "complete.obs"),
            p_value = cor_test(data = cur_data(), Ectomycorrhizal, AM, method = "pearson")$p)

cor_results <- guild_abundance %>%
  group_by(Fire.Severity) %>%
  summarise(correlation = cor(Ectomycorrhizal, AM, method = "pearson", use = "complete.obs"),
            p_value = cor_test(data = cur_data(), Ectomycorrhizal, Saprotroph, method = "pearson")$p)


print(cor_results)

library(lme4)
library(DHARMa)

hist(log10(guild_abundance$Ectomycorrhizal))
model_Interval<- lmer(log10(Ectomycorrhizal+1) ~ AM * Fire.Interval +(1/Site|Transect), data = guild_abundance)

model_Severity<- lmer(log10(Ectomycorrhizal+1) ~ AM * Fire.Severity +(1/Site|Transect), data = guild_abundance)

model_Regime<- lmer(log10(Ectomycorrhizal+1) ~ AM * (Fire.Interval + Fire.Severity) +(1/Site|Transect), data = guild_abundance)



sim_res <- simulateResiduals(fittedModel = model_Interval)
plot(sim_res)
#plot(explo_model_fire)
qqPlot(resid(model_Interval))
summary(model_Interval)
Anova_resin<-round(Anova(model_Interval,test='F'), 2) 
Anova_resin

sim_res <- simulateResiduals(fittedModel = model_Severity)
plot(sim_res)
#plot(explo_model_fire)
qqPlot(resid(model_Severity))
summary(model_Severity)
Anova_resin<-round(Anova(model_Severity,test='F'), 2) 
Anova_resin

sim_res <- simulateResiduals(fittedModel = model_Regime)
plot(sim_res)
#plot(explo_model_fire)
qqPlot(resid(model_Regime))
summary(model_Regime)
Anova_resin<-round(Anova(model_Regime,test='F'), 2) 
Anova_resin




ggplot(guild_abundance, aes(x = Saprotroph, y = AM, color = Fire.Interval)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Inter-guild interactions: ECM vs. Saprotroph",
       x = "Saprotrophic Fungi Abundance",
       y = "Ectomycorrhizal Fungi Abundance") +
  theme_minimal()

