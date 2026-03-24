



library(tidyverse)
library(ggplot2)
library(readxl)



myc_data <- read_excel("~/ABS_FIRE/ABS_FIRE_MYCO/Bag_data.xlsx", sheet=1)

myc_data$Tube_myc<-as.numeric(myc_data$Tube_myc)
myc_data$myc<-as.numeric(myc_data$myc)
myc_data$beads<-as.factor(myc_data$beads)
myc_data$Tube_ID<-as.factor(myc_data$Tube_ID)


bag_data <- read_excel("~/ABS_FIRE/ABS_FIRE_MYCO/Bag_data.xlsx", sheet=2)


bag_data$Location<-as.factor(bag_data$Location)
bag_data$Tube_ID<-as.factor(bag_data$Tube_ID)
bag_data$Site<-as.factor(bag_data$Site)
bag_data$Transect<-as.factor(bag_data$Transect)
bag_data$Nutrient_sub<-as.numeric(bag_data$Nutrient_sub)
bag_data$Bead_weight<-as.numeric(bag_data$Bead_weight)

bag_myc<-left_join(bag_data,myc_data, by='Tube_ID')%>%
  group_by(Site, Transect, Location) %>%
  mutate(Tube_ID = first(Tube_ID))%>%
  ungroup()

str(bag_data)
colnames(bag_data)


#compare myc data to myc with tube data

myc_data$myc_dif<-myc_data$Tube_myc-myc_data$Tube_w_mg

myc_data<-subset(myc_data,!Tube_ID%in% c('85','86'))

library(ggpmisc)
library(plotly)

#graph comaring just biomass vs biomass in tubes

p<-myc_data%>%
  ggplot(aes(x=myc,y=myc_dif))+
  geom_point(aes(color=beads),size=3)+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, size = 5) +
  labs(x = "Just_hyphae", y = "hyphae in tube") +
  theme_minimal()

ggplotly(p)
