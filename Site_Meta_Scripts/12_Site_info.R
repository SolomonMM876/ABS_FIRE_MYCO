library(tidyverse)


Bag_data<-read.csv('Processed_data/All_Bag_data.csv')

Bag_Sites<-Bag_data%>%
  select(Site,Fire.Severity,Fire.Interval, Vegetation.Class,Total.P,Carbon,Nitrogen,Ammonia_mg_kg,
         Nitrate_mg_kg,Ortho_P_mg_kg,pH,Days_Installed)%>%
  group_by(Site,Fire.Severity,Fire.Interval, Vegetation.Class)%>%
  summarise(round(across(Total.P:Days_Installed, ~ mean(.x, na.rm = TRUE)),digits=2), .groups = "drop")

head(Bag_Sites)
write.csv(Bag_Sites,"Processed_data/12_Sites_Info.csv", row.names=FALSE)
