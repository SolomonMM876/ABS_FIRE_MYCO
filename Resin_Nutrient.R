library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(ggplot2)
library(readxl)
library(tidyr)

Bag_Site_notes<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')%>%
  select(Site:mass_loss,-harvest_date)


#Ammonia and Nitrate data
NH4NO3_RESIN_LOT_1 <- read_csv("Raw_data/Nutes/24-07-10 SOLOMON NH4NO3 HCL RESIN LOT 1_edit.txt")
NH4NO3_RESIN_LOT_1$`Sample Details`[82]<-"10-2-3" #incorrectly typed in, should be site 10, not 1
# NH4NO3_RESIN_LOT_1<-NH4NO3_RESIN_LOT_1%>%
#   mutate(Result = ifelse(`Test Name` == "Nitrate 2" & Result < 0, 0.0757/2, Result)) #replace negative values with lowest possible/2


NH4NO3_RESIN_LOT_2 <- read_csv("Raw_data/Nutes/24-07-12 SOLOMON NH4NO3 RESIN HCL LOT 2_edit.txt")
NH4NO3_RESIN_LOT_2$`Sample Details`[17:18]<-"7-1-3" #only sample unaccounted for change unknown to ID
NH4NO3_RESIN_LOT_2$`Sample ID`[17:18]<-"7-1-3" #only sample unaccounted for 
# NH4NO3_RESIN_LOT_2<-NH4NO3_RESIN_LOT_2%>%
#   mutate(Result = ifelse(`Test Name` == "Nitrate 2" & Result < 0, 0.0587/2, Result)) #replace negative values with lowest possible/2

#lowest observable value 
LOQ_NH4NO3<-(.0757+.0587)/2

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
  select(-Sample_ID,-Sample_Details,-Time) %>% #
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
    # Blank_avg_Ammon = mean(Ammonia[grepl("CTRL", Site)], na.rm = TRUE),
    # Ammon_blanked = (Ammonia-Blank_avg_Ammon),
    Ammonia= ifelse(Ammonia < 0, LOQ_NH4NO3/2, Ammonia),#replace negative values with lowest possible/2
    Ammonia_mg_kg= Ammonia *(7.5/Nutrient_sub),
    #Nitrate_blanked = (Nitrate-Blank_avg_Nitrate),
    Nitrate = ifelse(Nitrate < 0, LOQ_NH4NO3/2, Nitrate), #replace negative values with lowest possible/2
    Nitrate_mg_kg= Nitrate*(7.5/Nutrient_sub))#7.5mL used for nutrient subset

#rm(NH4NO3_RESIN_LOT_1,NH4NO3_RESIN_LOT_2)

str(NO3_NH4_Resin)
NO3_NH4_Resin$Site = factor(NO3_NH4_Resin$Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56"))
NO3_NH4_Resin$Location = factor(NO3_NH4_Resin$Location, levels = c("3", "16","16.1","16.2", "33", "47"))


p<-NO3_NH4_Resin%>%
  filter(!grepl("CTRL", Site), na.rm = TRUE)%>%
  arrange(Site,Transect,Location)%>%
  ggplot(aes(x=Site, y=Ammonia_mg_kg ))+
  geom_col(aes(fill=Transect, color= Location), position = "dodge") +
  # facet_grid(~Fire.Interval,  scales = "free_x", )+
  #theme(legend.position = 'none')+
  labs(x = 'Site', y = 'Ammonia (mg/kg)') 

plotly::ggplotly(p)
p

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

# library(lubridate)
# NO3_NH4_Resin%>%
#   mutate(parsed_date = parse_date_time(`Time_Nitrate 2`, "a b d HMS Y"),
#          parsed_date = as.POSIXct(parsed_date))%>%
#    filter(!str_detect(Site, 'STAND|C C|CC|NIRAJ'))%>%
#   ggplot(aes(x=parsed_date, y=Nitrate ))+
#   geom_col(aes(group=Site,fill=Transect, color= Location), position = "dodge")







#Ortho-phos
OPHOS_RESIN_HCL_LOT_1 <- read_csv("Raw_data/Nutes/24-07-13 SOLOMON OPHOS RESIN HCL LOT 1.csv")
OPHOS_RESIN_HCL_LOT_2 <- read_csv("Raw_data/Nutes/24-07-13 SOLOMON OPHOS RESIN HCL LOT 2.csv")
OPHOS_RESIN_HCL_LOT_2$`Sample Details`[6]<-"29-1-3" #this sample needed to be manually adjusted- date issues



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
         Ortho_blanked = (Result-Blank_avg),
         Ortho_blanked = ifelse(Ortho_blanked < 0, 0.0288/2, Ortho_blanked),
         Ortho_P_mg_kg= Ortho_blanked *(7.5/Nutrient_sub)/(`Manual Dil`))#7.5mL used for 1 g of resin



Ortho_P_Resin$Site = factor(Ortho_P_Resin$Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56"))
Ortho_P_Resin$Location = factor(Ortho_P_Resin$Location, levels = c("3", "16","16.1","16.2", "33", "47"))


p<-Ortho_P_Resin%>%
  filter(!grepl("CTRL", Site), na.rm = TRUE)%>%
  arrange(Site,Transect,Location)%>%
  ggplot(aes(x=Site, y=Ortho_P_mg_kg ))+
  geom_col(aes(fill=Transect, color= Location), position = "dodge") +
 # facet_grid(~Fire.Interval,  scales = "free_x", )+
  #theme(legend.position = 'none')+
  labs(x = 'Site', y = 'PO4 (mg/kg)') 

plotly::ggplotly(p)

temp<-Ortho_P_Resin%>%
  group_by(Site,Transect)%>%
  summarize(n())

Nutrient_Resins<-inner_join(NO3_NH4_Resin%>%
                filter(!is.na(Site))%>%
                select(Site,Transect,Location, Ammonia_mg_kg, Nitrate_mg_kg,Tube_ID),
                Ortho_P_Resin%>%
                filter(!is.na(Site))%>%
                select(Site,Transect,Location, Ortho_P_mg_kg,Tube_ID))

library(writexl)
write_xlsx(Nutrient_Resins, path='Processed_data/Resin_Nutrients.xlsx')
            
