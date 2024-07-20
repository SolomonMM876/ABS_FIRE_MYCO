library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(ggplot2)

Bag_Site_notes<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')%>%
  select(Site:mass_loss,-harvest_date)





#Ammonia and Nitrate data
NH4NO3_RESIN_LOT_1 <- read_csv("Raw_data/Nutes/24-07-10 SOLOMON NH4NO3 HCL RESIN LOT 1.csv")
NH4NO3_RESIN_LOT_1$`Sample Details`[82]<-"10-2-3" #incorrectly typed in, should be site 10, not 1
NH4NO3_RESIN_LOT_2 <- read_csv("Raw_data/Nutes/24-07-12 SOLOMON NH4NO3 RESIN HCL LOT 2.csv")
NH4NO3_RESIN_LOT_2$`Sample Details`[17:18]<-"7-1-3" #only sample unaccounted for change unknown to ID
NH4NO3_RESIN_LOT_2$`Sample ID`[17:18]<-"7-1-3" #only sample unaccounted for 



NO3_NH4_Resin<-bind_rows(NH4NO3_RESIN_LOT_1,NH4NO3_RESIN_LOT_2)[,-c(10:13)] %>%
  rename(Test = `Test Name`, Sample_ID = `Sample ID`, Sample_Details = `Sample Details`) %>%
  #these case_when() select the appropriate Sample name
  mutate(Sample = case_when(
    !is.na(mdy(Sample_ID, quiet = TRUE)) ~ Sample_Details,
    TRUE ~ Sample_ID))%>%
  mutate(Sample = case_when(
    !is.na(dmy(Sample, quiet = TRUE)) ~ Sample_Details,
    TRUE ~ Sample))%>%
    mutate(Sample = case_when(
      is.na(Sample) ~ Sample_ID,
      TRUE ~ Sample))%>%
  select(-Time,-Sample_ID,-Sample_Details) %>%
  #this removes samples that are not needed 
  filter(!str_detect(Sample, 'STAND|C C|CC|NIRAJ'))%>%
  #this removes samples that are over diluted 
  filter(!str_detect(Sample, 'CTRL C|CTRL D'))%>%
  pivot_wider(names_from = Test, values_from = c(Result,Absorbance,`Auto Dil`)) %>%
  rename(Ammonia= 'Result_Ammonia 2.0',
         Nitrate = 'Result_Nitrate 2',
         Ammonia_Dil = `Auto Dil_Ammonia 2.0`)%>%
  separate(Sample, into = c("Site", "Transect", "Location"), sep = "-", convert = TRUE)%>%
  mutate(Transect=as.factor(Transect),
         Location=as.factor(Location))%>%
  #this mutate corrects small differences in site location so samples can join better
  mutate(Location = if_else(Location %in% c('03','2','5'), '3', Location),
         Location = if_else(Location %in% c('15','17'),'16',Location),
         Location = if_else(Location %in% c('32'),'33',Location))%>%
  left_join(Bag_Site_notes)%>%
  mutate(#calc blank avg for Nitr
    Blank_avg_Nitrate= mean(Nitrate[grepl("CTRL", Site)], na.rm = TRUE),
    Blank_avg_Ammon = mean(Ammonia[grepl("CTRL", Site)], na.rm = TRUE),
         Ammonia_mg_kg= (Ammonia-Blank_avg_Ammon)*(7.5/Nutrient_sub)/Ammonia_Dil,
         Nitrate_mg_kg= (Nitrate-Blank_avg_Nitrate)*(7.5/Nutrient_sub)/`Auto Dil_Nitrate 2`)#7.5mL used for nutrient subset

rm(NH4NO3_RESIN_LOT_1,NH4NO3_RESIN_LOT_2)


NO3_NH4_Resin$Site = factor(NO3_NH4_Resin$Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56"))
NO3_NH4_Resin$Location = factor(NO3_NH4_Resin$Location, levels = c("3", "16", "33", "47"))


p<-NO3_NH4_Resin%>%
  filter(!grepl("CTRL", Site), na.rm = TRUE)%>%
  arrange(Site,Transect,Location)%>%
  ggplot(aes(x=Site, y=Ammonia_mg_kg ))+
  geom_col(aes(fill=Transect, color= Location), position = "dodge") +
  # facet_grid(~Fire.Interval,  scales = "free_x", )+
  #theme(legend.position = 'none')+
  labs(x = 'Site', y = 'Ammonia (mg/kg)') 

plotly::ggplotly(p)

p<-NO3_NH4_Resin%>%
  filter(!grepl("CTRL", Site), na.rm = TRUE)%>%
  arrange(Site,Transect,Location)%>%
  ggplot(aes(x=Site, y=Nitrate_mg_kg ))+
  geom_col(aes(fill=Transect, color= Location), position = "dodge") +
  # facet_grid(~Fire.Interval,  scales = "free_x", )+
  #theme(legend.position = 'none')+
  labs(x = 'Site', y = 'Nitrate (mg/kg)') 

plotly::ggplotly(p)

temp<-NO3_NH4_Resin%>%
  group_by(Site,Transect)%>%
  summarize(n())




#Ortho-phos
OPHOS_RESIN_HCL_LOT_1 <- read_csv("Raw_data/Nutes/24-07-13 SOLOMON OPHOS RESIN HCL LOT 1.csv")
OPHOS_RESIN_HCL_LOT_2 <- read_csv("Raw_data/Nutes/24-07-13 SOLOMON OPHOS RESIN HCL LOT 2.csv")
OPHOS_RESIN_HCL_LOT_2$`Sample Details`[6]<-"29-1-3" #this sample needed to be manually adjusted- date issues
OPHOS_RESIN_HCL_LOT_2$`Sample Details`[10]<-"49-2-47" # manually adjusted extra 29-2-47 as this fit better
OPHOS_RESIN_HCL_LOT_2$`Sample ID`[10]<-"49-2-47" #
OPHOS_RESIN_HCL_LOT_2$`Sample Details`[9]<-"7-1-3" #only sample unaccounted for 
OPHOS_RESIN_HCL_LOT_2$`Sample ID`[9]<-"7-1-3" #only sample unaccounted for 



Ortho_P_Resin<-bind_rows(OPHOS_RESIN_HCL_LOT_1,OPHOS_RESIN_HCL_LOT_2) %>%
  rename(Test = `Test Name`, Sample_ID = `Sample ID`, Sample = `Sample Details`) %>%
  select(-Time,-Sample_ID, -Test)%>%
  filter(!str_detect(Sample, 'STAND|C C|CC|NIRAJ'))%>%
  separate(Sample, into = c("Site", "Transect", "Location"), sep = "-", convert = TRUE)%>%
  mutate(Transect=as.factor(Transect),
         Location=as.factor(Location))%>%
  mutate(Location = if_else(Location %in% c('03','2','5'), '3', Location),
         Location = if_else(Location %in% c('15','17'),'16',Location),
         Location = if_else(Location %in% c('32'),'33',Location))%>%
  filter(!grepl('REP',Location))%>% #these are reps of the same sample with simiilar values
  left_join(Bag_Site_notes)%>%
  mutate(Blank_avg = mean(Result[grepl("CTRL", Site)], na.rm = TRUE), 
         Ortho_P_mg_kg= (Result-Blank_avg)*(7.5/Nutrient_sub)/(`Manual Dil`*`Auto Dil`))#7.5mL used for 1 g of resin



Ortho_P_Resin$Site = factor(Ortho_P_Resin$Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56"))
Ortho_P_Resin$Location = factor(Ortho_P_Resin$Location, levels = c("3", "16", "33", "47"))


p<-Ortho_P_Resin%>%
  filter(!grepl("CTRL", Site), na.rm = TRUE)%>%
  arrange(Site,Transect,Location)%>%
  ggplot(aes(x=Site, y=Result ))+
  geom_col(aes(fill=Transect, color= Location), position = "dodge") +
 # facet_grid(~Fire.Interval,  scales = "free_x", )+
  #theme(legend.position = 'none')+
  labs(x = 'Site', y = 'PO4 (mg/kg)') 

plotly::ggplotly(p)

temp<-Ortho_P_Resin%>%
  group_by(Site,Transect)%>%
  summarize(n())

            