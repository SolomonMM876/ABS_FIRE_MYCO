library(lme4)
library(performance)
library(visreg)
#library(plotly)
library(DHARMa)
library(car)
library(emmeans)
library(readxl)
library(dplyr)

Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')
#both rows recorded 0 biomass because of harvest issue
#Removed row 14 because both bags were found out of the ground 
#Removed row 26 because no recorded biomass
x <- which(rowSums(is.na(select(Bag_Site, myc))) > 0)
Bag_Site[x, ]
Bag_Site <- Bag_Site[-x, ]  # drop those rows
Stoich_Totals <- read_excel("Raw_data/Stoich/Stoich_Totals_Round_1.xlsx")
#Nutrient resins
Resin_Nutrients<-read_excel('Processed_data/Resin_Nutrients.xlsx')

Bag_Site_Short<-Bag_Site%>%
  mutate(Regime = paste(Fire.Interval, Fire.Severity, sep = "_"))%>%
  dplyr::select(Site,Transect,Bray.P:Dead.Tree.Canopy.Cover_perc,Pair,log10_Second_Weight_bag_yield_est:Regime)%>%
  left_join(Stoich_Totals%>%rename(Carb_Hyph = Carbon, Nitrog_Hyph =Nitrogen, Phos_Hyph= Percent_Phos_))

Bag_data<-left_join(Bag_Site,Resin_Nutrients)%>%
  left_join(Stoich_Totals%>%rename(Carb_Hyph = Carbon, Nitrog_Hyph =Nitrogen, Phos_Hyph= Percent_Phos_))



#Resins related to biomass

# 
# m_biomass_resins<-lmer(log10_Second_Weight_bag_yield_est~  Ortho_P_mg_kg+ Nitrate_mg_kg+ 
#                          Fire.Interval + Fire.Severity+ (1|Site/Transect) , 
#                        data=Bag_data)
m_biomass_day_resins<-lmer(log10_biomass_day_all_cor~  Ortho_P_mg_kg+ Nitrate_mg_kg+ 
                         Fire.Interval + Fire.Severity+ (1|Site/Transect) , 
                       data=Bag_data)


summary(m_biomass_day_resins)
Anova_resin<-round(Anova(m_biomass_day_resins,test='F'), 2) 
Anova_resin
plot(m_biomass_day_resins)
qqPlot(resid(m_biomass_day_resins))
r2(m_biomass_day_resins)
emm_biomass_resins<-as.data.frame(emmeans(m_biomass_day_resins,
                                          ~Fire.Interval))


#Ortho P and Biomass
library(ggeffects)
predict_response(m_biomass_day_resins,terms= c('Ortho_P_mg_kg'), back_transform = FALSE)%>%
  mutate(Biomass_day= (10^(predicted)*(1e+06/15)),
         confidence.low= (10^(conf.low)*(1e+06/15)),
         confidence.high = (10^(conf.high)*(1e+06/15)))%>%
ggplot(aes(x, Biomass_day)) +
  geom_line(color = "black", linewidth=2) +
  geom_ribbon(aes(ymin =confidence.low , ymax = confidence.high), alpha = 0.1)+
  geom_point(data=Bag_data, mapping=aes(x=Ortho_P_mg_kg, y=(10^(log10_biomass_day_all_cor)*(1e+06/15))),
             inherit_aes=FALSE, size=3)+
  labs(
    x = expression(paste(PO[4], " (mg/kg)")),
    y = "Hyphal Production (g/ha/day)") +
  annotate("text", x = .5, y = Inf, label = paste0("Avail Phos (p) = ", Anova_resin$`Pr(>F)`[1]),
           hjust = .1, vjust = 1.1, size = 12)+
  theme_minimal(base_size = 15) +  # Minimal theme with larger base text size
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')



min(Bag_data$Biomass_day_all_cor)*(1e+06/15)

#Fire interval and biomass
interval_colors <- c("Long" = "darkred", "Short" = "orange")



p<-ggplot(emm_biomass_resins, aes(x = Fire.Interval, y = (10^(emmean)*(1e+06/15))) )+
  geom_col(aes(fill=Fire.Interval),size=4, width = .7) +
  geom_errorbar(aes(ymin = ((10^lower.CL)*(1e+06/15)),
                    ymax = ((10^upper.CL)*(1e+06/15)), width = 0.2 ),
                width= .4, size= 1.5) +
  geom_point(data=Bag_data, aes(x=Fire.Interval, y=(Biomass_day_all_cor)*(1e+06/15)), size=3, alpha=.6)+ 
  labs(x = "Fire Interval", y = "Hyphal Production (g/ha/day)") +
  scale_fill_manual(values = interval_colors) +     # Custom colors for Interval
  #scale_y_continuous(breaks = seq(0, 600, by = 100)) +
  annotate("text", x = 1.9, y = Inf, label = paste0("Interval (p) = ", Anova_resin$`Pr(>F)`[3]),
           hjust = 2.5, vjust = 1.5, size = 12)+
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

p



plotly::ggplotly(p)

emm_biomass_resins%>%
  summarise(emmeans_avg= mean(emmean))%>%
  mutate(emmeans_biomass_g_ha_day= (10^emmeans_avg)*(1e+06/15))

# in Bag_data_New.R script I calc the biomass production per day and the log of that on lines 164-172
Bag_Site%>%
  summarise(all_mean= mean(biomass_g_ha_day, na.rm=TRUE))


Bag_Site%>%
  ggplot(aes(x=Site,y=biomass_g_ha_day))+
  geom_col(aes(fill=Transect),position = "dodge")+
  facet_grid(~Fire.Severity,scales = "free_x")

severity_colors <- c("High" = "darkolivegreen", "Low" = "cornflowerblue")


p<-ggplot(emm_biomass_resins, aes(x = Fire.Severity, y = (10^(emmean)*(1e+06/15))) )+
  geom_col(aes(fill=Fire.Severity),size=4, width = .7) +
  geom_errorbar(aes(ymin = ((10^lower.CL)*(1e+06/15)),
                    ymax = ((10^upper.CL)*(1e+06/15)), width = 0.2 ),
                width= .4, size= 1.5) +
  labs(x = "Fire Severity", y = "Hyphal Production (g/ha/day)") +
  scale_fill_manual(values = severity_colors) +     # Custom colors for Interval
  scale_y_continuous(breaks = seq(0, 600, by = 100)) +
  annotate("text", x = 1.9, y = Inf, label = paste0("Severity (p) = ", Anova_resin$`Pr(>F)`[4]),
           hjust = 2.5, vjust = 1.5, size = 12)+
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')


p



#I trust the second round of weighing more, though both rounds produce similar results
m1<-lmer(log10_Second_Weight_bag_yield_est~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_Site)

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


#now adjusting for time that bags were in the ground
m1_day<-lmer(log10_biomass_day~  Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_Site)


#1
summary(m1_day)
Anova(m1_day,test='F')
plot(m1_day)
qqPlot(resid(m1_day))
pairs(emmeans(m1_day, ~ Fire.Interval))
emm_Interval<-as.data.frame(emmeans(m1_day, ~Fire.Interval))
emmip(m1_day, ~Fire.Interval)


#g/hectare:
#10,000 m^2 ×0.1 m= 1,000 m^3 
#1,000m^3 x (x mg /15cm^3) x (1g/1000mg) x 1000kg/ton x 1000g/kg= g/ha
(1e+06/15)


interval_colors <- c("Long" = "darkred", "Short" = "orange")



ggplot(emm_Interval, aes(x = Fire.Interval, y = (10^(emmean)*(1e+06/15))) )+
  geom_col(aes(fill=Fire.Interval),size=4, width = .7) +
  geom_errorbar(aes(ymin = ((10^lower.CL)*(1e+06/15)),
                    ymax = ((10^upper.CL)*(1e+06/15)), width = 0.2 ),
                width= .4, size= 1.5) +
  labs(x = "Fire Interval", y = "Hyphal Production (g/ha/day)") +
  scale_fill_manual(values = interval_colors) +     # Custom colors for Interval
  annotate("text", x = 1.5, y = Inf, label = "P = 0.083", hjust = 4, vjust = 1.1, size = 12)+
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')


r2(m1_day)

#2
m2<-lmer(log10(Ortho_P_mg_kg) ~   Fire.Interval + Fire.Severity + (1|Site) , data=Bag_data)
m3<-lmer(NH4~   Fire.Interval + Fire.Severity + (1|Site) , data=Bag_Site)
m4<-lmer(NO3~   Fire.Interval  + Fire.Severity+ (1|Site) , data=Bag_Site)
m5<-lmer(NO3~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)
m6<-lmer(NH4~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)
m7<-lmer(Bray.P~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)

hist(log10(Bag_Site_Short$C_N))
hist(sqrt(Bag_Site_Short$C_N), main = "Square Root Transformation", xlab = "sqrt(C_N)")
hist((Bag_Site_Short$C_N)^(1/3), main = "Cube Root Transformation", xlab = "Cube Root of C_N")
hist(1 / Bag_Site_Short$C_N, main = "Reciprocal Transformation", xlab = "1/C_N")# the best
hist(exp(Bag_Site_Short$C_N), main = "Exponential Transformation", xlab = "exp(C_N)")
z_score_transformed <- scale(Bag_Site_Short$C_N)
hist(z_score_transformed, main = "Z-score Normalization", xlab = "Z-score of C_N")


C_N_model<-lmer((1/C_N)~  Ortho_P_mg_kg+ Nitrate_mg_kg+ log10_biomass_day_all_cor +
                  Fire.Interval + Fire.Severity+ (1|Site) , data=Bag_data)
#8 C_N~ Fire Interval
summary(C_N_model)
Anova(C_N_model,test='F')
plot(C_N_model)
qqPlot(resid(C_N_model))
r2(C_N_model)
emm_C_N_Interval<-as.data.frame(emmeans(C_N_model, ~Fire.Interval))



interval_colors <- c("Long" = "darkred", "Short" = "orange")



p<-ggplot(emm_C_N_Interval, aes(x = Fire.Interval, y = 1/emmean) )+
  geom_col(aes(fill=Fire.Interval),size=4, width = .7) +
  geom_errorbar(aes(ymin = 1/lower.CL,
                    ymax = 1/upper.CL, width = 0.2 ),
                width= .4, size= 1.5) +
  labs(x = "Fire Interval", y = "Hyphal Carbon:Nitrogen") +
  scale_fill_manual(values = interval_colors) +     # Custom colors for Interval
 # annotate("text", x = 1.5, y = Inf, label = , hjust = 4, vjust = 1.1, size = 12)+
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')


plotly::ggplotly(p)

hist(log10(Bag_data$Ortho_P_mg_kg))



#2 (BrayP~Regime) Not sig
summary(m2)
Anova(m2,test='F')
plot(m2)
qqPlot(resid(m2))
r2(m2)
#3 (NH4~Regime) Not Sig
summary(m3)
Anova(m3,test='F')
plot(m3)
qqPlot(resid(m3))
pairs(emmeans(m3, ~ Fire.Severity+ Fire.Interval))
r2(m3)
#4 (NO3~Regime) 
#Severity almost sig t= -2.117, p= .06334
#Interval trend t= -1.807, p= .1042
#High Long - Low Short   df=8.02   t.ratio=2.824  p= 0.0850
summary(m4)
NO3_Anova<-Anova(m4,test='F')
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


grouped_Bag_Site <- Bag_Site %>%
  group_by(Fire.Interval) %>%
  summarise(
    mean_yield = mean(10^log10_Second_Weight_bag_yield_est, na.rm = TRUE),
    sd_yield = sd(10^log10_Second_Weight_bag_yield_est, na.rm = TRUE),
    count = n()
  )
###
dat_text <- data.frame(
  Fire.Interval = c('Long', 'Short'))

Bag_Site$Site = factor(Bag_Site$Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56"))


Bag_Site%>%
  ggplot(aes(x=Site, y=biomass_g_ha_day))+
  geom_col(aes(fill=Site), position = "dodge") +
  facet_grid(~Fire.Interval,  scales = "free_x", )+
theme(legend.position = 'none')+
    labs(x = 'Fire Interval', y = 'Biomass Production (mg)') ->p

plotly::ggplotly(p)
p

hist((Bag_Site$log10_Second_Weight_bag_yield_est))

Bag_Site %>%
  ggplot(aes(y = 10^(log10_biomass_day_all_cor),
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

Bag_Site %>%
  ggplot(aes(y = 10^(log10_Second_Weight_bag_yield_est), x = Site,)) +
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

Bag_Site %>%
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




Bag_Site %>%
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
pairs(~ log10_biomass_day_all_cor+Nitrate_mg_kg + Ortho_P_mg_kg + Herb.Cover_0.50cm_perc:Live.Tree.Canopy.Cover_perc  , data = Bag_data)






# Load the library for stepwise selection
library(MASS)

# Full model with all potential explanatory variables
full_model <- lm(myc_bag_yield_est ~ ., data = na.omit(Bag_Site[,c(24,27:34,36:44,46:62)]))

# Perform stepwise selection
stepwise_model <- stepAIC(full_model, direction = "both", trace = FALSE)
summary(stepwise_model)

# Fit the null model and other candidate models
null_model <- lm(Corrected_myc ~ 1, data = Bag_Site)

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

