
library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)

Bag_Site<-read_excel('~/ABS_FIRE/ABS_FIRE_MYCO/Processed_data/All_Bag_Site_Info.xlsx')


Myc_Weight<-Bag_Site%>%
  select(Site,Transect,Location,myc)%>%
  mutate(Location_Group = case_when(
    Location %in% c(3, 16) ~ "3_and_16",
    Location %in% c(33, 47) ~ "33_and_47"))%>%
  group_by(Site,Transect, Location_Group)%>%
    mutate(sum_myc = sum(myc),
           whats_left= sum_myc-2.5)


Myc_Weight%>%
  ggplot(aes(x=Site, y=sum_myc))+
  geom_col(aes(fill=Transect, color= Location_Group), position = "dodge")+
  geom_text(aes(label=sum_myc))


#One subsample (~1–3 mg) was analysed for hyphal C and N concentration


#For P the other subsample (~1–2 mg)

96/12
