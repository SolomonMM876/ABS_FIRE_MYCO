library(tidyverse)
library(lme4)
library(performance)
library(visreg)
library(plotly)
library(DHARMa)
library(car)
library(emmeans)
library(readxl)

Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')
Myc_Weight<-Bag_Site%>%
  select(Site,Transect,Location,myc)%>%
  mutate(Location_Group = case_when(
    Location %in% c(3, 16) ~ "3_and_16",
    Location %in% c(33, 47) ~ "33_and_47"))%>%
  group_by(Site,Transect, Location_Group)%>%
  mutate(sum_myc = sum(myc),
         whats_left= sum_myc-2.5)

#both rows recorded 0 biomass because of harvest issue
#Removed row 14 because both bags were found out of the ground
#Removed row 26 because no recorded biomass
#or do below and drop whole row
x <- which(rowSums(select(Bag_Site, myc)) == 0)
Bag_Site[x, ]
Bag_Site <- Bag_Site[-x, ]  # drop those rows


# here's where I started editing code
Bag_Site %>%
  # sum biomass by transect end
  mutate(Location_Group = case_when(
    Location %in% c(3, 16, 16.1, 16.2) ~ "3_and_16",
    Location %in% c(33, 47) ~ "33_and_47"))%>%
  group_by(Site,Transect, Location_Group)%>%
  mutate(sum_myc = sum(myc),
         sum_by_lg = sum(10^log10_myc_bag_yield_est), 
         whats_left= sum_myc-2.5, 
         ) %>% 
  ungroup() %>% 
  # sum biomass by groups (lowest and highest biomass samples)
  group_by(Site,Transect) %>% 
  arrange(log10_myc_bag_yield_est, .by_group=TRUE) %>% 
  mutate(rank = 1:n(), 
         Amount_Group = case_when(rank %in% 1:2 ~ 'low', 
                                  rank %in% 3:4 ~ 'high')) %>% 
  group_by(Site, Transect, Amount_Group) %>% 
  mutate(sum_by_amt = sum(10^log10_myc_bag_yield_est)) %>% 
  # sum biomass by transect
  group_by(Site,Transect) %>% 
  mutate(sum_by_transect = sum(10^log10_myc_bag_yield_est)) -> Bag_Site1  

# area plot to look at mass along each transect
Bag_Site1 %>% 
  mutate(Location = as.numeric(Location)) %>% 
  ggplot(aes(x=Location, y=10^log10_myc_bag_yield_est)) +
  geom_area() +
  facet_grid(rows=vars(Site), cols=vars(Transect))

# original model
m1<-lmer(log10_myc_bag_yield_est~Fire.Severity+ Fire.Interval + (1|Site/Transect), data=Bag_Site1)
summary(m1)
Anova(m1, test='F')

# model pooling by transect ends
tmp <- Bag_Site1 %>% 
  group_by(Site, Transect, Location_Group) %>% 
  slice(1) %>% 
  ungroup()
m2<-lmer(log10(sum_by_lg)~Fire.Severity+ Fire.Interval + (1|Site/Transect), data=tmp)
summary(m2)
Anova(m2, test='F')

# model pooling by sample mass
tmp <- Bag_Site1 %>% 
  group_by(Site, Transect, Amount_Group) %>% 
  slice(1) %>% 
  ungroup()
m3<-lmer(log10(sum_by_amt)~Fire.Severity+ Fire.Interval + (1|Site/Transect), data=tmp)
summary(m3)
Anova(m3, test='F')

# model pooling by transect
tmp <- Bag_Site1 %>% 
  group_by(Site, Transect) %>% 
  slice(1) %>% 
  ungroup()
m4<-lmer(log10(sum_by_amt)~Fire.Severity+ Fire.Interval + (1|Site), data=tmp)
summary(m4)
Anova(m4, test='F')
