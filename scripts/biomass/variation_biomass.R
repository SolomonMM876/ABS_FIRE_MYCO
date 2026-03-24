library(dplyr)
library(ggplot2)




Bag_Site<-read_excel('~/ABS_FIRE/ABS_FIRE_MYCO/processed_data/All_Bag_Site_Info.xlsx')
Bag_Site$Fire.Interval<-factor(Bag_Site$Fire.Interval, levels = c('Short','Long'))
Bag_Site$Fire.Severity<-factor(Bag_Site$Fire.Severity, levels = c('Low','High'))


Bag_Site$log10_myc_bag_yield_est <- log10(Bag_Site$myc_bag_yield_est + 0.18)


hist(log10(.05+Bag_Site$NO3))

#removing outliers
Bag_Site$log10_myc_bag_yield_est[c(14,26)]<-NA
Bag_Site$myc_bag_yield_est[c(14,26)]<-NA
#both rows recorded 0 biomass because of harvest issue
#Removed row 14 because both bags were found out of the ground 
#Removed row 26 because no recorded biomass


variation_by_site<-Bag_Site%>%
  group_by(Site,Fire.Interval) %>%
  summarise(mean_yield = mean(log10_myc_bag_yield_est,na.rm = TRUE), sd_yield = sd(log10_myc_bag_yield_est,na.rm = TRUE)) %>%
  mutate(cov_yield = sd_yield/mean_yield)%>%
  arrange(desc(cov_yield))%>%
  ungroup()

t_test_result <- t.test(cov_yield ~ Fire.Interval, data = variation_by_site)
print(t_test_result)



# Perform F-test to compare variances
f_test_result <- var.test(Bag_Site$log10_myc_bag_yield_est[Bag_Site$Fire.Interval == "Long"], 
                          Bag_Site$log10_myc_bag_yield_est[Bag_Site$Fire.Interval == "Short"])
print(f_test_result)


# Perform Levene's test to compare variances
levene_test_result <- leveneTest(log10_myc_bag_yield_est ~ Fire.Interval, data = Bag_Site)
print(levene_test_result)

# Boxplot for visualization
ggplot(Bag_Site, aes(x = Fire.Interval, y = log10_myc_bag_yield_est, fill = Fire.Interval)) +
  geom_boxplot() +
  labs(title = "Yield by Fire Interval", x = "Fire Interval", y = "Yield") +
  theme_minimal()








m_cov<-lmer(cov_yield~ Fire.Interval + (1|Site) , data=variation_by_site)

summary(m_cov)
Anova(m_cov,test='F')
plot(m_cov)
qqPlot(resid(m_cov))
r2(m_cov)



