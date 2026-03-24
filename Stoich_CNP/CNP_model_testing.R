library(lme4)
library(tidyr)
library(DHARMa)
library(car)
library(performance)
library(emmeans)
library(ggplot2)

#samples 11.1 and 56.2 are weird, try analyses with and without
#Sample 56.2 had lower than normal values for CNH, suggesting there was something in the sample besides hyphae
#one possibility could be resin beads that were not removed

#run DF joins first
source('df_joins.R')

Bag_data<-read.csv('Processed_data/All_Bag_data.csv')

CNP_clean<-read.csv("Processed_data/CNP_clean.csv")

Bag_data<-Bag_data%>%
  group_by(Site, Transect) %>%
  summarise(Ortho_P_mg_kg = mean(Ortho_P_mg_kg, na.rm = TRUE),
            Nitrate_mg_kg = mean(Nitrate_mg_kg, na.rm = TRUE),
            Ammonia_mg_kg = mean(Ammonia_mg_kg, na.rm = TRUE)) %>% 
  left_join(CNP_clean)



#CN#######
C_N_model<-lmer(log10(C_N)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                  Fire.Interval + Fire.Severity+ 
                  (1|Site) , 
                data = Bag_data )# %>% filter(!Sample %in% c("56.2")))
summary(C_N_model)
Anova(C_N_model,test='F')
plot(C_N_model)
qqPlot(resid(C_N_model))
#r2(C_N_model)
emm_C_N_Ortho<-as.data.frame(emmeans(C_N_model, ~Ortho_P_mg_kg))
emm_C_N_Ortho
emm_C_N_Severity<-as.data.frame(emmeans(C_N_model, ~Fire.Severity))
emm_C_N_Severity



interval_colors <- c("Long" = "darkred", "Short" = "orange")

library(ggeffects)
predict_response(C_N_model,terms= c('Fire.Severity'))
  
emm_C_N_Severity%>%  
ggplot(aes(x = Fire.Severity, y = 10^emmean) )+
  geom_boxplot(aes(fill=Fire.Severity),linewidth=4, width = .7) +
  geom_errorbar(aes(ymin = 10^lower.CL,
                    ymax = 10^upper.CL, width = 0.2 ),
                width= .4, size= 1.5) +
  labs(x = "Fire Severity", y = "Hyphal Carbon:Nitrogen") +
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')

Bag_data%>%
  ggplot(aes(x=as.factor(Site)) )+ 
  geom_point(aes(y=C_N), size = 5, shape= 'square') +
  geom_boxplot(aes(x=Fire.Severity, y= C_N))+
  facet_grid(~Fire.Severity, scales = 'free_x')+
  labs( y= ('C:N'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

slope <- fixef(C_N_model)["Ortho_P_mg_kg"]
intercept <- fixef(C_N_model)["(Intercept)"]
# Plot
Bag_data %>%
  ggplot(aes(x = Ortho_P_mg_kg, y = C_N)) + 
  geom_point(size = 5, shape = 'circle') +  # Raw data points
  labs(y = 'C:N', x = 'Ortho Phos') +
  geom_abline(slope = slope, intercept = intercept, color = "black", linewidth = 2) +  # Regression line
  geom_ribbon(data = emm_C_N_Ortho,inherit.aes = FALSE, aes(x = Ortho_P_mg_kg, ymin = lower.CL, ymax = upper.CL), alpha = 0.1) +  # Confidence interval
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

Bag_data


library(glmmTMB)
#CP model########


hist(Bag_data$C_P)
hist(log10(Bag_data$C_P))


C_P_model<-lmer(log10(C_P)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                  Fire.Interval + Fire.Severity+
                  (1|Site), 
                  data=Bag_data)

sim_res <- simulateResiduals(fittedModel = C_P_model)
plot(sim_res)
# Test for overdispersion
testDispersion(sim_res)

# Test for uniformity of residuals
testUniformity(sim_res)


summary(C_P_model)
Anova(C_P_model)
plot(C_P_model)
qqPlot(resid(C_P_model))
#r2(C_P_model)
emm_C_N_Interval<-as.data.frame(emmeans(C_P_model, ~Fire.Interval))
emm_C_N_Severity<-as.data.frame(emmeans(C_P_model, ~Fire.Severity))


#NP model
N_P_model<-lmer(log10(N_P)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                      Fire.Interval + Fire.Severity+ 
                  (1|Site) , data=Bag_data)

sim_res <- simulateResiduals(fittedModel = N_P_model)
plot(sim_res)
#8 C_N~ Fire Interval
summary(N_P_model)
Anova(N_P_model,test='F')
plot(N_P_model)
qqPlot(resid(N_P_model))
#r2(N_P_model)
emm_C_N_Interval<-as.data.frame(emmeans(N_P_model, ~Fire.Interval))
emm_C_N_Severity<-as.data.frame(emmeans(N_P_model, ~Fire.Severity))



#single elements
Carb_Hyph_model<-lmer((Carbon)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                  Fire.Interval + Fire.Severity+ 
                  (1|Site) , 
                data = Bag_data)
summary(Carb_Hyph_model)
Anova(Carb_Hyph_model,test='F')
plot(Carb_Hyph_model)
qqPlot(resid(Carb_Hyph_model))

Nitrog_Hyph_model<-lmer((Nitrogen)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                        Fire.Interval + Fire.Severity+ 
                        (1|Site) , 
                      data = Bag_data)# %>% filter(!Sample %in% c("56.2")))
summary(Nitrog_Hyph_model)
Anova(Nitrog_Hyph_model,test='F')
plot(Nitrog_Hyph_model)
qqPlot(resid(Nitrog_Hyph_model))


Phos_Hyph_model<-lmer((Percent_Phos_ )~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                          Fire.Interval + Fire.Severity+ 
                          (1|Site) , 
                      data=Bag_data %>% filter(!Sample %in% c("11.1","12.2")))

summary(Phos_Hyph_model)
Anova(Phos_Hyph_model,test='F')
plot(Phos_Hyph_model)
qqPlot(resid(Phos_Hyph_model))

