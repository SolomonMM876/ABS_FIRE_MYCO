library(lme4)
library(performance)
library(visreg)
library(plotly)
library(DHARMa)
library(car)
library(emmeans)
library(readxl)
library(dplyr)

Bag_Site<-read_excel('~/ABS_FIRE/ABS_FIRE_MYCO/Processed_data/All_Bag_Site_Info.xlsx')


Tubes<-Bag_Site%>%
  select(Site,Transect,Location,Tube_ID,myc)%>%
  group_by(Site,Transect)
library(writexl)
write_xlsx(Tubes,'~/ABS_FIRE/ABS_FIRE_MYCO/Processed_data/Tube_ID.xlsx')

Myc_Weight<-Bag_Site%>%
  select(Site,Transect,Location,myc)%>%
  mutate(Location_Group = case_when(
    Location %in% c(3, 16) ~ "3_and_16",
    Location %in% c(33, 47) ~ "33_and_47"))%>%
  group_by(Site,Transect, Location_Group)%>%
  mutate(sum_myc = sum(myc),
         whats_left= sum_myc-2.5)




#removing outliers
Bag_Site$log10_myc_bag_yield_est[c(14,26)]<-NA
Bag_Site$myc_bag_yield_est[c(14,26)]<-NA
#both rows recorded 0 biomass because of harvest issue
#Removed row 14 because both bags were found out of the ground 
#Removed row 26 because no recorded biomass
#or do below and drop whole row
x <- which(rowSums(select(Bag_Site, myc)) == 0)
Bag_Site[x, ]
Bag_Site <- Bag_Site[-x, ]  # drop those rows



m1<-lmer(log10_myc_bag_yield_est~Fire.Severity+ Fire.Interval + (1|Site/Transect) , data=Bag_Site)
m2<-lmer(Bray.P~   Fire.Interval + Fire.Severity + (1|Site) , data=Bag_Site)
m3<-lmer(NH4~   Fire.Interval + Fire.Severity + (1|Site) , data=Bag_Site)
m4<-lmer(NO3~   Fire.Interval  * Fire.Severity+ (1|Site) , data=Bag_Site)
m5<-lmer(NO3~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)
m6<-lmer(NH4~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)
m7<-lmer(Bray.P~   log10_myc_bag_yield_est + (1|Site) , data=Bag_Site)




#1
summary(m1)
#this is how much more the short fire interval increases myco biomass
(10^0.15274)-.05
Anova(m1,test='F')
#residual vs fitted plot 
plot(m1)
qqPlot(resid(m1))
# Bag_Site_Outliers<-Bag_Site[c(14,26),]
pairs(emmeans(m1, ~ Fire.Severity+ Fire.Interval))
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
Anova(m4,test='F')
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


Bag_Site%>%
  ggplot(aes(x=Site, y=myc_bag_yield_est))+
  geom_col(aes(fill=Transect, color= Location), position = "dodge") +
  facet_grid(~Fire.Interval,  scales = "free_x", )+
  # geom_text(aes(label = beads, y = myc_bag_yield_est + 0.05, group = interaction(Location,Transect)), position = position_dodge(width = 1), 
  #           color = 'black',  size = 4)+
theme(legend.position = 'none')
#alpha = mass_loss
#+
  # # Add mean points
  # geom_point(data = variation_by_site, aes(x = as.factor(Site), y = mean_yield), 
  #            color = 'red', size = 3) +
  # # Add error bars for standard deviation
  # geom_errorbar(data = variation_by_site, 
  #               aes(x = as.factor(Site),y=mean_yield, ymin = mean_yield - sd_yield, 
  #                   ymax = mean_yield + sd_yield), 
  #               width = 0.2, color = 'red')
ggplotly(p)

hist(log10(Bag_Site$myc_bag_yield_est))

Bag_Site%>%
  ggplot(aes(x=Fire.Interval, y= myc_bag_yield_est))+
  geom_col(aes(Fill=Fire.Interval))+theme_minimal()






# Pairwise scatter plots to see relationships
pairs(~ myc_bag_yield_est + FESM_severity +Bray.P + NH4 + NO3 + Total.P + Carbon + Nitrogen, data = Bag_Site)






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

