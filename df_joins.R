library(tidyr)
library(readxl)
library(dplyr)

Bag_Site<-read.csv('Processed_data/All_Bag_Site_Info.csv')
#both rows recorded 0 biomass because of harvest issue
#Removed row 14 because both bags were found out of the ground 
#Removed row 26 because no recorded biomass
x <- which(rowSums(is.na(select(Bag_Site, myc_2nd_w ))) > 0)
Bag_Site[x, ]
Bag_Site <- Bag_Site[-x, ]  # drop those rows


Stoich_Totals <- read.csv("Processed_data/CNP_clean.csv")%>%
  mutate(Site=as.numeric(Site),
         Transect=as.numeric(Transect))
#Nutrient resins
Resin_Nutrients<-read.csv('Processed_data/Soil_pH_w_resin_nutes.csv')#soil pH is also here

Myco_abun<-read.csv('Processed_data/Myco_host_abundance.csv')
N_fixer_abun<-read.csv('Processed_data/N_fixer_abundance.csv')
#seeder_resprout_abun<-read.csv('Processed_data/resprout_seeder_abundance.csv')


#biomass data corrected for seq
#Biomass_seq_corrected<-read_excel('Processed_data/corrected_biomass_df.xlsx')

Bag_data<-left_join(Bag_Site,Resin_Nutrients)%>%
  mutate(Regime = paste(Fire.Interval, Fire.Severity, sep = "_"))%>%
  left_join(Stoich_Totals%>%rename(Carb_Hyph = Carbon, Nitrog_Hyph =Nitrogen, Phos_Hyph= Percent_Phos_))%>%
  #calc the average time between fire 
  mutate(
    diff1 = Most.Recent.Fire_Year - `X2nd.Most.Recent.Fire_Year`,
    diff2 = `X2nd.Most.Recent.Fire_Year` - as.numeric(`X3rd.Most.Recent.Fire_Year`),
    fire.frequency = rowMeans(cbind(diff1, diff2), na.rm = TRUE))%>%
  left_join(Myco_abun)%>%
  left_join(N_fixer_abun)
  #%>%left_join(fire_response_abun)



write.csv(Bag_data, file='Processed_data/All_Bag_data.csv',row.names=FALSE)

