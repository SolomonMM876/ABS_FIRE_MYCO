library(lme4)
library(performance)
library(visreg)
library(plotly)
library(DHARMa)
library(car)
library(emmeans)
library(readxl)
library(dplyr)

Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')

Bag_Site_Short<-Bag_Site%>%
  mutate(Regime = paste(Fire.Interval, Fire.Severity, sep = "_"))%>%
  dplyr::select(Site,Transect,Bray.P:Dead.Tree.Canopy.Cover_perc,Pair,log10_Second_Weight_bag_yield_est:Regime)




  #removing outliers
  #Bag_Site$log10_myc_bag_yield_est[c(14,26)]<-NA
  #Bag_Site$myc_bag_yield_est[c(14,26)]<-NA
#both rows recorded 0 biomass because of harvest issue
#Removed row 14 because both bags were found out of the ground 
#Removed row 26 because no recorded biomass
#or do below and drop whole row
x <- which(rowSums(select(Bag_Site, myc)) == 0)
Bag_Site[x, ]
Bag_Site <- Bag_Site[-x, ]  # drop those rows



m1<-lmer(log10_Second_Weight_bag_yield_est~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_Site)
#checking to see if second round of weighing changed anything
m2<-lmer(Bray.P~   Fire.Interval + Fire.Severity + (1|Site) , data=Bag_Site)
m3<-lmer(NH4~   Fire.Interval + Fire.Severity + (1|Site) , data=Bag_Site)
m4<-lmer(NO3~   Fire.Interval  + Fire.Severity+ (1|Site) , data=Bag_Site)
m5<-lmer(NO3~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)
m6<-lmer(NH4~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)
m7<-lmer(Bray.P~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)




#1
summary(m1)
#this is how much more the short fire interval increases myco biomass
(10^0.13700)-0.095
Anova(m1,test='F')
#residual vs fitted plot 
plot(m1)
qqPlot(resid(m1))
# Bag_Site_Outliers<-Bag_Site[c(14,26),]
pairs(emmeans(m1, ~ Fire.Interval))
emm_Interval<-as.data.frame(emmeans(m1, ~Fire.Interval))
emmip(m1, ~Fire.Interval)
ggplot(emm_Interval, aes(x = Fire.Interval, y = 10^emmean+.095)) +
  geom_point(aes(),size=4) +
  geom_errorbar(aes(ymin = 10^lower.CL+.095, ymax = 10^upper.CL+.095), width = 0.2 ) +
  labs(x = "Fire Interval", y = "Estimated Marginal Mean", title = "Estimated Marginal Means by Fire Interval") +
    annotate("text", x = 1.5, y = Inf, label = "P = 0.066", hjust = 1.1, vjust = 1.1, size = 5)+
  theme_minimal()

r2(m1)
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
#8 (Bray.P~Biomass)
summary(m8)
Anova(m8,test='F')
plot(m8)
qqPlot(resid(m8))
r2(m8)

grouped_Bag_Site <- Bag_Site %>%
  group_by(Fire.Interval) %>%
  summarise(
    mean_yield = mean(10^log10_Second_Weight_bag_yield_est, na.rm = TRUE),
    sd_yield = sd(10^log10_Second_Weight_bag_yield_est, na.rm = TRUE),
    count = n()
  )
###
dat_text <- data.frame(
  label = c("Fire Interval = Pr(>F) 0.066 ", ""),
  Fire.Interval = c('Long', 'Short'))

Bag_Site$Site = factor(Bag_Site$Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56"))


Bag_Site%>%
  ggplot(aes(x=Site, y=10^log10_Second_Weight_bag_yield_est ))+
  geom_col(aes(fill=Transect, color= Location), position = "dodge") +
  facet_grid(~Fire.Interval,  scales = "free_x", )+
theme(legend.position = 'none')+
    labs(x = 'Fire Interval', y = 'Biomass Production (mg)') 

# # Add mean points
# geom_point(data = variation_by_site, aes(x = as.factor(Site), y = mean_yield),
#            color = 'red', size = 3) +
# # Add error bars for standard deviation
# geom_errorbar(data = variation_by_site,
#               aes(x = as.factor(Site),y=mean_yield, ymin = mean_yield - sd_yield,
#                   ymax = mean_yield + sd_yield),
#               width = 0.2, color = 'red')
ggplotly(p)

hist((Bag_Site$log10_Second_Weight_bag_yield_est))

Bag_Site %>%
  ggplot(aes(y = 10^log10_Second_Weight_bag_yield_est+.095, x = Fire.Interval)) +
  geom_boxplot(aes(), width = 0.1, fill = c('Dark Red', 'Orange')) +
  geom_point(aes(fill = Fire.Interval), alpha = 0.5, size = 2) +
  annotate("text", x = 0.5, y = Inf, label = 'Fire Interval = Pr(>F) 0.066 ', 
           vjust = 1.5, hjust = -3, size = 4) + # Adjust x, y, vjust, and hjust to position the text
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
  ggplot(aes(y = 10^log10_Second_Weight_bag_yield_est+.095, x = Site,)) +
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
pairs(~ myc_bag_yield_est+Bray.P + NH4 + NO3 + Total.P + Carbon + Nitrogen, data = Bag_Site)






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

