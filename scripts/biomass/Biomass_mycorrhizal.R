library(lme4)
library(performance)
library(visreg)
library(ggplot2)
library(DHARMa)
library(car)
library(emmeans)
library(readxl)
library(tidyverse)
library(tibble)

Bag_data<-read.csv('Processed_data/All_Bag_data.csv')



#read in modeling_bags_2nd_rnd

source("ABS_Second_Rnd/modeling_bags_2nd_rnd.R")
#From Bag_ITS_Prep
guild_summary<- read.csv('Processed_data/Percent_myco_biomass')

temp<-Bag_data%>%
  left_join(guild_summary)

ggplot(temp, aes(x =ratio_myco_nonmyco  , y =biomass_g_ha_day )) +
  geom_point(color = "blue", alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    x = "Mycorrhizal:Total Reads",
    y = 'Biomass g ha day '  ) +
  theme_minimal()

temp_filter<-temp %>% 
  filter(ratio_myco_nonmyco>.5) %>% 
  filter(!is.na(log10_biomass_day))

fitler<-temp %>% 
  filter(ratio_myco_nonmyco<.5) 

temp_adj<-temp %>% 
  mutate(myco_biomass=biomass_g_ha_day*ratio_myco_nonmyco,
         log_myco_biomass = ifelse(myco_biomass == 0, NA, log10(myco_biomass)))

#Resins related to biomass
# 
# m_biomass_day_resins<-lmer(log_myco_biomass~Fire.Severity+Fire.Interval +
#                              Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg + #soil prop
#                              perc_myco_host_freq+ 
#                              (1|Site/Transect) 
#                            ,data=temp_adj)


m_additive <- lmer(log10_biomass_day ~ Fire.Severity + Fire.Interval +
                     Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + pH +
                     perc_myco_host_freq +
                     (1 | Site/Transect), data = temp_filter)


m_interaction <- lmer(log10_biomass_day ~ Fire.Severity * Fire.Interval +
                        Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + pH +
                        perc_myco_host_freq +
                        (1 | Site/Transect), data = temp_filter)


AIC(m_additive, m_interaction)



m_biomass_day_resins<-lmer(log10_biomass_day~Fire.Severity+Fire.Interval +
                             Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg + pH+ #soil prop
                             perc_myco_host_freq+
                             (1|Site/Transect) 
                           ,data=temp_filter)

#################################
#making a table of relationships between factors and bioamss

library(broom.mixed)

# Tidy fixed effect estimates
fixed_effects <- broom.mixed::tidy(m_biomass_day_resins, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Predictor = gsub("Low|Short", "", term)               # Clean labels
  )

# Type II ANOVA table
anova_df <- as.data.frame(Anova(m_biomass_day_resins, test = "F"))
anova_df$Predictor <- rownames(anova_df)

# Join slope estimates with p-values
summary_table <- fixed_effects %>%
  rename(Slope = estimate, Std_Error = std.error) %>%
  left_join(anova_df[, c("Predictor", "Pr(>F)",'Df', 'Df.res','F')], by = "Predictor") %>%
  rename(p = `Pr(>F)`) %>% 
  mutate(    Predictor = gsub("_mg_kg", "", Predictor),                  # Optional: simplify nutrient names
             Predictor = gsub("_", " ", Predictor),                      # Optional: make labels more readable
             Predictor = trimws(Predictor))

# Add empty R2 column to model summary table, then bind
summary_table_1st <- summary_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  select(-effect,-term,-statistic) %>% 
  mutate(Round='1st') %>% 
  relocate(Round,Predictor,F,Slope,Std_Error,Df,Df.res,p) 

#Bind results from first and second round together

summary_table<-rbind(summary_table_1st,summary_table_2nd)

library(flextable)
library(officer)

# Build flextable
ft_biomass <- summary_table %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Linear mixed-effects model of mycorrhizal biomass by fire factors and soil nutrients. Significant terms in bold (p ≤ 0.05).") %>%
  bold(i = ~ p <= 0.05, j = c("Slope", "p"), bold = TRUE) %>%
  hline(i = nrow(summary_table), border = fp_border(width = 1), part = "body") %>%
  hline(i = nrow(summary_table_1st), border = fp_border(width = 1), part = "body") %>%
  align(align = "left", part = "all") %>%
  autofit()

# Show flextable in viewer
ft_biomass

# Save to Word document
doc <- read_docx() %>%
  body_add_flextable(ft_biomass)

print(doc, target = "Tables/Biomass_Model_Summary_Table.docx")








############################################
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


#MAKE a table#####
temp<-Anova_resin%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>%
  mutate(Round='First Round')

temp
write.csv(temp,'Tables/Biomass_Anova.csv', row.names = FALSE)

mean(temp_filter$biomass_g_ha_day)*(150/1000)

emm_Interval <- emm_biomass_Interval %>%
  mutate(
    emmean_bt = 10^emmean*(1e+06/30)* (150/1000),
    lower.CL_bt = 10^lower.CL*(1e+06/30)* (150/1000),
    upper.CL_bt = 10^upper.CL*(1e+06/30)* (150/1000)
  )
emm_Interval

interval_colors <- c("Long" = "darkred", "Short" = "orange")


p<-ggplot(emm_biomass_Interval, aes(x = Fire.Interval, y = (10^(emmean)*(1e+06/30)* (150/1000))) )+
  geom_boxplot(aes(fill=Fire.Interval),size=4, width = .7) +
  geom_jitter(data=Bag_data, aes(x=Fire.Interval, y=(Biomass_day)*(1e+06/30)* (150/1000)), size=3,width = .05, alpha=.6)+
  geom_errorbar(aes(ymin = ((10^lower.CL)*(1e+06/30)* (150/1000)),
                    ymax = ((10^upper.CL)*(1e+06/30)* (150/1000)) ),
                width= .2, size= 1.5) +
  labs(x = "Fire Frequency", y = "Hyphal Production (kg/ha) over 5 months") +
  scale_fill_manual(values = interval_colors) +     # Custom colors for Interval
  scale_y_continuous(breaks = seq(0, 350, by = 50)) +
  annotate("text", x =1.75, y = Inf, label = paste0("Fire Frequency (p) = ", Anova_resin["Fire.Interval", "Pr(>F)"]),
           hjust = 2.5, vjust = 1.5, size = 10)+
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(linewidth = 2, colour = "black"),
        legend.position = 'none')
p

# Back-transform the EMMs
emm_Ortho <- emm_biomass_Ortho %>%
  arrange(Ortho_P_mg_kg) %>%  # sort by x-axis
  mutate(
    emmean_bt = 10^emmean * (1e+06 / 30) * (150 / 1000),
    lower.CL_bt = 10^lower.CL * (1e+06 / 30) * (150 / 1000),
    upper.CL_bt = 10^upper.CL * (1e+06 / 30) * (150 / 1000)
  )
emm_Ortho
# Plot
#Ortho P and Biomass
library(ggeffects)
predict_response(m_biomass_day_resins,terms= c('Ortho_P_mg_kg'), back_transform = FALSE)%>%
  mutate(Biomass_day= (10^(predicted)*(1e+06/30)* (150 / 1000)),
         confidence.low= (10^(conf.low)*(1e+06/30)* (150 / 1000)),
         confidence.high = (10^(conf.high)*(1e+06/30)* (150 / 1000)))%>%
  ggplot(aes(x, Biomass_day)) +
  geom_line(color = "black", linewidth=2) +
  geom_ribbon(aes(ymin =confidence.low , ymax = confidence.high), alpha = 0.1)+
  geom_point(data=temp_filter, mapping=aes(x=Ortho_P_mg_kg, y=(10^(log10_biomass_day)*(1e+06/30)* (150 / 1000))),
     size=3)+
  labs(
    x = expression(paste(PO[4], " (mg/kg)")),
    y = "Hyphal Production (kg/ha/day)") +
  annotate("text", x = 0, y = Inf, label = paste0("Avail Phos (p) = ", Anova_resin$`Pr(>F)`[3]),
           hjust = .1, vjust = 1.1, size = 12)+
  theme_minimal(base_size = 15) +  # Minimal theme with larger base text size
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')->p

plotly::ggplotly(p)
####################################
###########plots by sample location 
#2 different ways


#####biomass actually recovered

temp_adj %>%
  group_by(site_transect_locat = paste("Site", Site, "T", Transect, 'L',Location)) %>%
  mutate(
    Site = as.factor(Site),
    site_transect = factor(site_transect_locat,
                           levels = unique(site_transect_locat)[order(Site, Transect,Location)]))%>%
  ggplot(aes(x = site_transect_locat, y = myc_2nd_w)) +
  geom_bar(stat = "identity", color='black') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Site, scales = "free") +
  labs(
    x = "Site/Transect",
    y = "recovered biomass (ug)",
    fill = "Guild" ) 

#####biomass estimated recovered

temp_adj %>%
  group_by(site_transect_locat = paste("Site", Site, "T", Transect, 'L',Location)) %>%
  mutate(
    Site = as.factor(Site),
    site_transect = factor(site_transect_locat,
                           levels = unique(site_transect_locat)[order(Site, Transect,Location)]))%>%
  ggplot(aes(x = site_transect_locat, y = myc_2nd_w_est_yield)) +
  geom_bar(stat = "identity", color='black') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Site, scales = "free") +
  labs(
    x = "Site/Transect",
    y = "adjusted biomass  (ug)",
    fill = "Guild" ) 
