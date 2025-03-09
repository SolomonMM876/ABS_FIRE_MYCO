library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(ggpubr)

Nutrients_Transects<-read_excel('Processed_data/Nutrients_Transect_level.xlsx')
Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')

Twelve_sites <- semi_join( Nutrients_Transects,Bag_Site, by = c("Transect", "Site"))%>%
  arrange(as.numeric(Site))%>%
  left_join(Bag_Site%>%select(Site,Transect,Fire.Interval))

Twelve_sites$Site<-as.factor(Twelve_sites$Site)


 Twelve_long<- Twelve_sites %>%
  pivot_longer(
    cols = Bray.P:Nitrogen,
    names_to = "Nutrient",
    values_to = "Value" )%>%
  arrange(as.numeric(Site))
 
facet_order <- c("NO3", "NH4", "Nitrogen", "Bray.P", "Total.P", "Carbon")
 
 
Twelve_long%>%
ggplot( aes(x = Site, y = Value)) +
  geom_boxplot() +
  facet_wrap(~ Nutrient, scales = "free_y") +
  labs(x = "Site", y = "Nutrient Value (mg/kg)") +
  theme_minimal()

# Arrange the data by Site
Twelve_sites <- Twelve_sites %>%
  mutate(Site = factor(Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56")))

# Create each plot individually
plot_NO3 <- ggplot(Twelve_sites, aes(x = Site, y = NO3)) +
  geom_boxplot() +
  labs(x = NULL, y = "NO3 (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_NH4 <- ggplot(Twelve_sites, aes(x = Site, y = NH4)) +
  geom_boxplot() +
  labs(x = NULL, y = "NH4 (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_Nitrogen <- ggplot(Twelve_sites, aes(x = Site, y = Nitrogen)) +
  geom_boxplot() +
  labs(x = NULL, y = "Nitrogen (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_BrayP <- ggplot(Twelve_sites, aes(x = Site, y = Bray.P)) +
  geom_boxplot() +
  labs(x = NULL, y = "Bray.P (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_TotalP <- ggplot(Twelve_sites, aes(x = Site, y = Total.P)) +
  geom_boxplot() +
  labs(x = "Site", y = "Total.P (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_Carbon <- ggplot(Twelve_sites, aes(x = Site, y = Carbon)) +
  geom_boxplot() +
  labs(x = "Site", y = "Carbon (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

# Arrange plots using ggarrange and remove specific elements
ggarrange(
  plot_NO3 + rremove("x.text") + rremove("x.axis"),
  plot_NH4 + rremove("x.text") + rremove("x.axis"),
  plot_Nitrogen + rremove("x.text") + rremove("x.axis"),
  plot_BrayP + rremove("x.text") + rremove("x.axis"),
  plot_TotalP + rremove("y.text") + rremove("x.top.axis"),
  plot_Carbon + rremove("y.text") + rremove("y.axis"),
  ncol = 2, nrow = 3
)
