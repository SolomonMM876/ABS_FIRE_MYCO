library(tidyr)
library(dplyr)
library(readxl)
library(stringr)

Ph_Soil<- read_excel("Raw_data/Ph_kumari_Dec_24.xlsx")
Nutrient_Resins<-read.csv('Processed_data/Resin_Nutrients.csv')%>%
  mutate(Transect=as.factor(Transect),
         Location=as.factor(Location))




Ph_Soil<-Ph_Soil%>%
  mutate(
    Site = str_extract(Sample_ID, "(?<=S)\\d+"),
    Transect = str_extract(Sample_ID, "(?<=T)\\d+"),
    Location = str_extract(Sample_ID, "(?<=L).+")
  )%>%
  mutate(Transect=as.factor(Transect),
         Location=as.factor(Location))%>%
  #this mutate corrects small differences in site location so samples can join better
  mutate(Location = if_else(Location %in% c('03','2','5'), '3', Location),
         Location = if_else(Location %in% c('15','17'),'16',Location),
         Location = if_else(Location %in% c('32'),'33',Location),
         Location = if_else(Location %in% c('48'),'47',Location))

Ph_Resin<-right_join(Ph_Soil%>%mutate(Site=as.factor(Site)),Nutrient_Resins%>%mutate(Site=as.factor(Site)))



Ph_Resin<-Ph_Resin%>%
  select(-Lab_No,-Sample_ID,-soil_wt_g)

write.csv(Ph_Resin, file='Processed_data/Soil_pH_w_resin_nutes.csv',row.names=FALSE)

