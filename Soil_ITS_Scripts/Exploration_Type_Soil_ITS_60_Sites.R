library(ggplot2)
library(lme4)
source('Soil_ITS_Scripts/Soil_ITS_data prep.R')
library(ggplot2)





dat_explo<-dat_myco_RA%>%
 #rownames_to_column("SH_ID_OTU")%>%
  filter(guild=='mycorrhizal')%>%
  group_by(Site,Transect,exploration_type)%>%
  summarise(explo_count=sum(readcount))%>%
  left_join(dat_myco_RA%>%select(Site,Transect,Severity,Interval,reads_samp)%>%distinct())%>%
  mutate(RA_explo=explo_count/reads_samp,
         log_RA_explo=log10(RA_explo+0.00028/2))
#this is relative abundance compared to other mycorrhizal fungi

    


hist(dat_explo$RA_explo)
hist(log10(dat_explo$RA_explo))

library(glmmTMB)
library(DHARMa)


explo_model_fire<-lmer((log_RA_explo)~ exploration_type *(Severity+Interval) + 
                            (1|Site/Transect),
                          data=dat_explo)


#fire model
sim_res <- simulateResiduals(fittedModel = explo_model_fire)
plot(sim_res)
#plot(explo_model_fire)
qqPlot(resid(explo_model_fire))
summary(explo_model_fire)
Anova_resin<-round(Anova(explo_model_fire,test='F'), 2) 
Anova_resin


#GRaphing sig results
emm_just_Interval <- emmeans(explo_model_fire, ~ Interval)

emm_severity <- emmeans(explo_model_fire, ~ exploration_type * Severity)
emm_interval <-emmeans(explo_model_fire, ~ exploration_type * Interval)
pairs_severity<-as.data.frame( pairs(emm_severity, adjust = "tukey"))  # Tukey HSD
pairs_interval<-as.data.frame( pairs(emm_interval, adjust = "tukey")) # Tukey HSD
emm_severity_df <-as.data.frame( emmeans(explo_model_fire, ~ exploration_type * Severity))
emm_interval_df <-as.data.frame( emmeans(explo_model_fire, ~ exploration_type * Interval))

#Interval
# Filter only significant contrasts (p ≤ 0.05)
sig_pairs <- pairs_interval %>%
  filter(p.value <= 0.05) %>%
  select(contrast, p.value)


# Extract the group names from the "contrast" column
sig_pairs <- sig_pairs %>%
  mutate(
    group1 = gsub(" ", ".", gsub("[()]", "", sub(" - .*", "", contrast))),  # Extract first group and replace space with '.'
    group2 = gsub(" ", ".", gsub("[()]", "", sub(".* - ", "", contrast)))   # Extract second group and replace space with '.'
  )

sig_pairs

filtered_sig_pairs <- sig_pairs %>%
  separate(group1, into = c("prefix1", "suffix1"), sep = "\\.") %>%
  separate(group2, into = c("prefix2", "suffix2"), sep = "\\.") %>%
  filter(prefix1 == prefix2, suffix1 != suffix2) %>%
  unite("group1", prefix1, suffix1, sep = ".") %>%
  unite("group2", prefix2, suffix2, sep = ".")

comparison_list <- as.list(as.data.frame(t(filtered_sig_pairs[, c("group1", "group2")])))


library(ggsignif)
p<-emm_interval_df%>%
  ggplot( aes(x = interaction(exploration_type, Interval), y = (10^(emmean))) )+
  geom_col(size=3, alpha=.6)+
  geom_errorbar(aes(ymin=(10^(emmean-SE)),ymax=(10^(emmean+SE))), width=.2) +
  labs(x = "Fire Interval + Exploration Type", y = "RA_exploration_type") +
  theme_classic()+
  geom_signif(
    comparisons = comparison_list,
    annotations = formatC(filtered_sig_pairs$p.value, digits = 2),
    y_position = seq(max(10^(emm_interval_df$emmean)) * 1.4, 
                     max(10^(emm_interval_df$emmean)) * 1.1, 
                     length.out = nrow(filtered_sig_pairs)),  # Stagger y positions  # Second comparison (higher) 
    tip_length = 0.03, 
    vjust = 0.8, 
    size = 1.2) +
  theme(axis.text.x = element_text( hjust = 1, size = 15, face = "bold",  angle = 45),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

p
########################
#Severity
# Filter only significant contrasts (p ≤ 0.05)
sig_pairs <- pairs_severity %>%
  filter(p.value <= 0.05) %>%
  select(contrast, p.value)

# Extract the group names from the "contrast" column
sig_pairs <- sig_pairs %>%
  mutate(
    group1 = gsub(" ", ".", gsub("[()]", "", sub(" - .*", "", contrast))),  # Extract first group and replace space with '.'
    group2 = gsub(" ", ".", gsub("[()]", "", sub(".* - ", "", contrast)))   # Extract second group and replace space with '.'
  )

filtered_sig_pairs <- sig_pairs %>%
  separate(group1, into = c("prefix1", "suffix1"), sep = "\\.") %>%
  separate(group2, into = c("prefix2", "suffix2"), sep = "\\.") %>%
  filter(prefix1 == prefix2, suffix1 != suffix2) %>%
  unite("group1", prefix1, suffix1, sep = ".") %>%
  unite("group2", prefix2, suffix2, sep = ".")

comparison_list <- as.list(as.data.frame(t(filtered_sig_pairs[, c("group1", "group2")])))


library(ggsignif)
p<-emm_severity_df%>%
  ggplot( aes(x = interaction(exploration_type, Severity), y = (10^(emmean))) )+
  geom_col(size=3, alpha=.6)+
  geom_errorbar(aes(ymin=(10^(emmean-SE)),ymax=(10^(emmean+SE))), width=.2) +
  labs(x = "Fire Severity + Exploration Type", y = "RA_exploration_type") +
  theme_classic()+
  geom_signif(
    comparisons = comparison_list,
    annotations = formatC(filtered_sig_pairs$p.value, digits = 2),
    y_position = seq(max(10^(emm_severity_df$emmean)) * 1.4, 
                     max(10^(emm_severity_df$emmean)) * 1.1, 
                     length.out = nrow(filtered_sig_pairs)),  # Stagger y positions  # Second comparison (higher) 
    tip_length = 0.03, 
    vjust = 0, 
    size = 1.2) +
  theme(axis.text.x = element_text( hjust = 1, size = 15, face = "bold",  angle = 45),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

p  

# #nutrient model
# sim_res <- simulateResiduals(fittedModel = explo_model_fire)
# plot(sim_res)
# #plot(explo_model_fire)
# qqPlot(resid(explo_model_fire))
# summary(explo_model_fire)
# Anova_resin<-round(Anova(explo_model_fire), 2) 
# Anova_resin

