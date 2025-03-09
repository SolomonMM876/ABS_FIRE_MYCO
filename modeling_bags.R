library(lme4)
library(performance)
library(visreg)
library(ggplot2)
library(DHARMa)
library(car)
library(emmeans)
library(readxl)
library(dplyr)

Bag_data<-read.csv('Processed_data/All_Bag_data.csv')



#log transform response variables
Bag_data<-Bag_data%>%
  mutate(Ortho_P_mg_kg_log=log10(Ortho_P_mg_kg),
         Nitrate_mg_kg_log=log10(Nitrate_mg_kg),
         Ammonia_mg_kg_log=log10(Ammonia_mg_kg),
         pH_log=log10(pH))



#Resins related to biomass

m_biomass_day_resins<-lmer(log10_biomass_day~Fire.Severity+Fire.Interval +
                                 Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg + pH+ #soil prop
                                 perc_myco_host_freq+
                                 (1|Site/Transect) 
                               ,data=Bag_data)


#anova(m_biomass_day_resins_freq,m_biomass_day_resins_tot)

summary(m_biomass_day_resins)
Anova_resin<-round(Anova(m_biomass_day_resins,test='F'), 2) 
Anova_resin
plot(m_biomass_day_resins)
qqPlot(resid(m_biomass_day_resins))
r2(m_biomass_day_resins)
emm_biomass_Interval<-as.data.frame(emmeans(m_biomass_day_resins,
                                          ~Fire.Interval))
emm_biomass_Ortho<-as.data.frame(emmeans(m_biomass_day_resins,
                                            ~Ortho_P_mg_kg))
emm_Interval<-as.data.frame(emmeans(m3, ~Fire.Interval))
emm_Interval <- emm_Interval %>%
  mutate(
    emmean_bt = 10^emmean,
    lower.CL_bt = 10^lower.CL,
    upper.CL_bt = 10^upper.CL
  )



#Resins related to biomass

m_biomass_day_nutri<-lmer(log10_biomass_day~Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg + pH+ #soil prop
                             (1|Site/Transect) 
                           ,data=Bag_data)


#anova(m_biomass_day_resins_freq,m_biomass_day_resins_tot)

summary(m_biomass_day_nutri)
Anova_resin<-round(Anova(m_biomass_day_nutri,test='F'), 2) 
Anova_resin
plot(m_biomass_day_nutri)
qqPlot(resid(m_biomass_day_nutri))
r2(m_biomass_day_nutri)
emm_biomass_Interval<-as.data.frame(emmeans(m_biomass_day_nutri,
                                            ~Fire.Interval))


#Ortho P and Biomass
library(ggeffects)
predict_response(emm_biomass_Ortho,terms= c('Ortho_P_mg_kg'), back_transform = FALSE)%>%
  mutate(Biomass_day= (10^(predicted)*(1e+06/15)),
         confidence.low= (10^(conf.low)*(1e+06/15)),
         confidence.high = (10^(conf.high)*(1e+06/15)))%>%
ggplot(aes(x, Biomass_day)) +
  geom_line(color = "black", linewidth=2) +
  geom_ribbon(aes(ymin =confidence.low , ymax = confidence.high), alpha = 0.1)+
  geom_point(data=Bag_data, mapping=aes(x=myco_host_total, y=(10^(lbiomass_day)*(1e+06/15))),
             inherit_aes=FALSE, size=3)+
  labs(
    x = expression(paste(PO[4], " (mg/kg)")),
    y = "Hyphal Production (g/ha/day)") +
  annotate("text", x = 0, y = Inf, label = paste0("Avail Phos (p) = ", Anova_resin$`Pr(>F)`[1]),
           hjust = .1, vjust = 1.1, size = 12)+
  theme_minimal(base_size = 15) +  # Minimal theme with larger base text size
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')



min(Bag_data$Biomass_day)*(1e+06/15)

Fire interval and biomass
emm_biomass_Interval%>%
  summarise(emmeans_avg= mean(emmean))%>%
  mutate(emmeans_biomass_g_ha_day= (10^emmeans_avg)*(1e+06/15))

# in Bag_data_New.R script I calc the biomass production per day and the log of that on lines 164-172
# Bag_data%>%
#   summarise(all_mean= mean(biomass_g_ha_day, na.rm=TRUE))
# 
# 
# Bag_data%>%
#   ggplot(aes(x=Site,y=biomass_g_ha_day))+
#   geom_col(aes(fill=Transect),position = "dodge")+
#   facet_grid(~Fire.Severity,scales = "free_x")

severity_colors <- c("High" = "darkolivegreen", "Low" = "cornflowerblue")
interval_colors <- c("Long" = "darkred", "Short" = "orange")


p<-ggplot(emm_biomass_Interval, aes(x = Fire.Interval, y = (10^(emmean)*(1e+06/15))) )+
  geom_col(aes(fill=Fire.Interval),size=4, width = .7) +
  geom_point(data=Bag_data, aes(x=Fire.Interval, y=(Biomass_day)*(1e+06/15)), size=3, alpha=.6)+
  geom_errorbar(aes(ymin = ((10^lower.CL)*(1e+06/15)),
                    ymax = ((10^upper.CL)*(1e+06/15)), width = 0.2 ),
                width= .4, size= 1.5) +
  labs(x = "Fire Frequency", y = "Hyphal Production (g/ha/day)") +
  scale_fill_manual(values = interval_colors) +     # Custom colors for Interval
  scale_y_continuous(breaks = seq(0, 2300, by = 250)) +
  annotate("text", x = 1.6, y = Inf, label = paste0("Severity (p) = ", Anova_resin["Fire.Interval", "Pr(>F)"]),
           hjust = 2.5, vjust = 1.5, size = 12)+
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')


p


#Orthophosphate
m1<-lmer(log10(Ortho_P_mg_kg)~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_data)

#m1<-lmer((Ortho_P_mg_kg_log)~Fire.Severity+ Fire.Interval + Nitrate_mg_kg+ Ammonia_mg_kg + pH+ (1|Site/Transect) , data=Bag_data)



#1
summary(m1)
Anova(m1,test='F')
plot(m1)
qqPlot(resid(m1))
pairs(emmeans(m1, ~ Fire.Interval))
emm_Interval<-as.data.frame(emmeans(m1, ~Fire.Interval))
emmip(m1, ~Fire.Interval)
r2(m1)

#g/hectare 10,000 m^2 ×0.1 m= 1,000 m^3 
# 1,000m^3 x (x mg /15cm^3) x (1g/1000mg) x 1000kg/ton= kg/ha
(1000/15)


#Nitrate
m2<-lmer(log10(Nitrate_mg_kg)~  Fire.Severity+ Fire.Interval +(1|Site/Transect) , data=Bag_data)


#m2<-lmer((Nitrate_mg_kg_log)~Fire.Severity+ Fire.Interval + Ortho_P_mg_kg + Ammonia_mg_kg + pH+ (1|Site/Transect) , data=Bag_data)


#1
summary(m2)
Anova(m2,test='F')
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
m3<-lmer(log10(Ammonia_mg_kg) ~   Fire.Interval + Fire.Severity + (1|Site/Transect) , data=Bag_data)


#m3<-lmer((Ammonia_mg_kg_log)~Fire.Severity+ Fire.Interval + Ortho_P_mg_kg + Nitrate_mg_kg + pH+ (1|Site/Transect) , data=Bag_data)

#3 (NH4~Regime) Not Sig
summary(m3)
Anova(m3,test='F')
plot(m3)
qqPlot(resid(m3))
r2(m3)
emm_Interval<-as.data.frame(emmeans(m3, ~Fire.Interval))
emm_Interval <- emm_Interval %>%
  mutate(
    emmean_bt = 10^emmean,
    lower.CL_bt = 10^lower.CL,
    upper.CL_bt = 10^upper.CL
  )


#host identity
m4<-lmer(perc_myco_host_freq~   Fire.Interval  + Fire.Severity+ (1|Site) , data=Bag_data)




m5<-lmer(NO3~   log10_myc_2nd_w_est_yield + (1|Site) , data=Bag_data)
m6<-lmer(NH4~   log10_myc_2nd_w_est_yield + (1|Site) , data=Bag_data)
m7<-lmer(Bray.P~   log10_myc_2nd_w_est_yield + (1|Site) , data=Bag_data)


hist(log10(Bag_data$Ortho_P_mg_kg))





#4 (Host~Regime

summary(m4)
Anova<-Anova(m4,test='F')
Anova
plot(m4)
qqPlot(resid(m4))
pairs(emmeans(m4, ~ Fire.Severity+ Fire.Interval))

r2(m4)
#5 (No3~Biomass)
summary(m5)
Anova(m5,test='F')
plot(m5)
qqPlot(resid(m5))
r2(m5)
#6 (NH4~Biomass)
summary(m6)
Anova(m6,test='F')
plot(m6)
qqPlot(resid(m6))
r2(m6)
#7 (Bray.P~Biomass)
summary(m7)
Anova(m7,test='F')
plot(m7)
qqPlot(resid(m7))
r2(m7)


grouped_Bag_data <- Bag_data %>%
  group_by(Fire.Interval) %>%
  summarise(
    mean_yield = mean(10^log10_myc_2nd_w_est_yield, na.rm = TRUE),
    sd_yield = sd(10^log10_myc_2nd_w_est_yield, na.rm = TRUE),
    count = n()
  )
###
dat_text <- data.frame(
  Fire.Interval = c('Long', 'Short'))

Bag_data$Site = factor(Bag_data$Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56"))


Bag_data%>%
  ggplot(aes(x=Site, y=biomass_g_ha_day))+
  geom_col(aes(fill=Site), position = "dodge") +
  facet_grid(~Fire.Interval,  scales = "free_x", )+
theme(legend.position = 'none')+
    labs(x = 'Fire Interval', y = 'Biomass Production (mg)') ->p

plotly::ggplotly(p)
p

hist((Bag_data$log10_Second_Weight_bag_yield_est))

Bag_data %>%
  ggplot(aes(y = 10^(log10_myc_2nd_w_est_yield),
             x = Fire.Interval)) +
  geom_boxplot(aes(), width = 0.1, fill = c('Dark Red', 'Orange')) +
  geom_point(aes(fill = Fire.Interval), alpha = 0.5, size = 2) +
 # annotate("text", x = 0.5, y = Inf, label = 'Fire Interval = Pr(>F) 0.066 ', 
       #    vjust = 1.5, hjust = -3, size = 4) + # Adjust x, y, vjust, and hjust to position the text
  labs(x = 'Fire Interval', y = 'Biomass Production (mg)') +
  theme( 
    axis.title.x = element_text(size = rel(1.5)),   
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.text.x = element_text(size = rel(2)), # Increase the size of x-axis text
    axis.text.y = element_text(size = rel(1.5))
  ) +
   guides(fill="none")+
  theme_minimal()

Bag_data %>%
  ggplot(aes(y = 10^(log10_myc_2nd_w_est_yield), x = Site,)) +
  geom_boxplot(aes(fill= Fire.Interval), width = 0.1) +
    scale_fill_manual(values = c("Long" = "darkred", "Short" = "orange")) +
  geom_point(aes(fill = Fire.Interval), alpha = 0.5, size = 2) +
  facet_grid(~Pair, scales = "free_x") +
  annotate("text", x = 0.5, y = Inf, label = 'Fire Interval = Pr(>F) 0.066 ', 
           vjust = 1.5, hjust = -3, size = 4) + # Adjust x, y, vjust, and hjust to position the text
  labs(x = 'Fire Interval', y = 'Biomass Production (mg)') +
  theme( 
    axis.title.x = element_text(size = rel(1.5)),   
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.text.x = element_text(size = rel(2)), # Increase the size of x-axis text
    axis.text.y = element_text(size = rel(1.5))
  ) +
   #guides(fill="none")+
  theme_minimal()

Bag_data %>%
  ggplot(aes(y = 10^log10_Second_Weight_bag_yield_est + 0.095, x = Site, fill = Fire.Interval)) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = c("Long" = "darkred", "Short" = "orange")) +
  geom_point(alpha = 0.5, size = 2) +
  facet_grid(~Pair, scales = "free_x") +
  annotate("text", x = 0.5, y = Inf, label = 'Fire Interval = Pr(>F) 0.066', 
           vjust = 1.5, hjust = -0.1, size = 4) +
  labs(x = 'Site', y = 'Biomass Production (mg)') +
  theme_minimal() +
  theme( 
    axis.title.x = element_text(size = rel(1.5)),   
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.text.x = element_text(size = rel(2)), 
    axis.text.y = element_text(size = rel(1.5))
  ) +
  geom_smooth(method = "lm", aes(group = Pair), se = FALSE, color = "blue") +
  stat_smooth(method = "lm", aes(group = Pair), se = FALSE, color = "blue", 
              geom = "text", label = after_stat(expr(paste("R^2 = ", round(summary(lm(y ~ x))$r.squared, 2)))), 
              vjust = -0.5, size = 3)




Bag_data %>%
  ggplot(aes(y = NO3, x =Fire.Interval))+
  geom_boxplot(aes(fill = Fire.Interval), width = 0.2) +
  geom_point(aes(fill = Fire.Interval), alpha = 0.5, size = 2) +
 scale_fill_manual(values = c('Long' = 'Dark Red', 'Short' = 'Orange')) +
  #facet_grid(~Regime, scales = "free_x") +
   # geom_text(data = NO3_Anova, mapping = aes(x = 0, y = Inf,label = `Pr(>F)`),
   #           vjust = 1.5, hjust = -0.1, size = 3) + # Specify the facet for annotation
  labs(x = 'Site', y = 'NO3 mg/kg') +
  theme(
    axis.title.x = element_text(size = rel(1.5)),   
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.text.x = element_text(size = rel(2)), 
    axis.text.y = element_text(size = rel(1.5)),
    legend.position = "top" # Remove the legend
  ) +
    guides(fill="none")+
  theme_minimal()



# Pairwise scatter plots to see relationships
pairs(~ log10_myc_2nd_w_est_yield+Nitrate_mg_kg + Ortho_P_mg_kg + Herb.Cover_0.50cm_perc:Live.Tree.Canopy.Cover_perc  , data = Bag_data)






# Load the library for stepwise selection
library(MASS)

# Full model with all potential explanatory variables
full_model <- lm(myc_bag_yield_est ~ ., data = na.omit(Bag_data[,c(24,27:34,36:44,46:62)]))

# Perform stepwise selection
stepwise_model <- stepAIC(full_model, direction = "both", trace = FALSE)
summary(stepwise_model)

# Fit the null model and other candidate models
null_model <- lm(Corrected_myc ~ 1, data = Bag_data)

# Compare AIC 
AIC(null_model, full_model, stepwise_model)

# Load library for VIF
library(car)

# Calculate VIF for the stepwise selected model
vif(full_model)




# Site.info<-read_excel('Site.Info.xlsx')
# Site.info$Transect<-as.factor(Site.info$Transect)
# Site.info<-Site.info%>%
#   mutate(Site= gsub('ABS00|ABS0', "", Site))

