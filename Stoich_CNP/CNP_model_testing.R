library(lme4)
library(tidyr)
library(DHARMa)
library(car)
library(performance)
library(emmeans)
library(ggplot2)

#samples 11.1 and 56.2 are weird, try analyses with and without

#run DF joins first
source('df_joins.R')

Bag_data<-read.csv('Processed_data/All_Bag_data.csv')



hist(log10(Bag_data$C_N))
hist(sqrt(Bag_data$C_N), main = "Square Root Transformation", xlab = "sqrt(C_N)")
hist((Bag_data$C_N)^(1/3), main = "Cube Root Transformation", xlab = "Cube Root of C_N")
hist(1 / Bag_data$C_N, main = "Reciprocal Transformation", xlab = "1/C_N")# the best
hist(exp(Bag_data$C_N), main = "Exponential Transformation", xlab = "exp(C_N)")

#CN#######
C_N_model<-lmer((C_N)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                  Fire.Interval + Fire.Severity+ 
                  (1|Site) , 
                data = Bag_data)# %>% filter(!Sample %in% c("56.2")))
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
ggplot(aes(x = Fire.Severity, y = emmean) )+
  geom_boxplot(aes(fill=Fire.Severity),linewidth=4, width = .7) +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL, width = 0.2 ),
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

library(glmmTMB)
#CP model########


hist(Bag_data$C_P)
hist(log10(Bag_data$C_P))


C_P_model<-lmer(log10(C_P)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                  Fire.Interval + Fire.Severity+ log10_biomass_day+
                  (1|Site), 
                  data=Bag_data %>% filter(!Sample %in% c("11.1")))

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
                  log10_biomass_day +
                  Fire.Interval + Fire.Severity+ 
                  (1|Site) , data=Bag_data%>% filter(!Sample %in% c('11.1',"56.2")))

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
Carb_Hyph_model<-lmer((Carb_Hyph)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                  Fire.Interval + Fire.Severity+ 
                  (1|Site) , 
                data = Bag_data)# %>% filter(!Sample %in% c("56.2")))
summary(Carb_Hyph_model)
Anova(Carb_Hyph_model,test='F')
plot(Carb_Hyph_model)
qqPlot(resid(Carb_Hyph_model))

Nitrog_Hyph_model<-lmer((Nitrog_Hyph)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                        Fire.Interval + Fire.Severity+ 
                        (1|Site) , 
                      data = Bag_data)# %>% filter(!Sample %in% c("56.2")))
summary(Nitrog_Hyph_model)
Anova(Nitrog_Hyph_model,test='F')
plot(Nitrog_Hyph_model)
qqPlot(resid(Nitrog_Hyph_model))


Phos_Hyph_model<-lmer((Phos_Hyph)~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
                          Fire.Interval + Fire.Severity+ 
                          (1|Site) , 
                        data = Bag_data)# %>% filter(!Sample %in% c("56.2")))
summary(Phos_Hyph_model)
Anova(Phos_Hyph_model,test='F')
plot(Phos_Hyph_model)
qqPlot(resid(Phos_Hyph_model))
