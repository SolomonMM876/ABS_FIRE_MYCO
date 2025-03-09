
#LOAD LIBRARARIES

library(tidyverse)
library(vegan)
library(ggplot2)




#All meta data from 12 sites with bags collected
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')[,-1]

Bag_Site_ID<-read.csv('Processed_data/All_Bag_Site_Info.csv')%>%
  filter(!Tube_ID %in% c(8,11))#these tubes are both from site 29, but I mixed up the transects and 
#it creates duplication errors when I keep samples in
#I could include them in Site level analyses.

Bag_data<-Bag_data%>%left_join(Bag_Site_ID, relationship = "many-to-many")
rm(Bag_Site_ID)


Seq_data<-read.csv('Processed_data/cleaned_seq_dat.csv')%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))%>%
  filter(guild %in%c("Ectomycorrhizal","Arbuscular Mycorrhizal","Orchid Mycorrhizal","Ericoid Mycorrhizal"))%>%
  select(-sample)

tax<-Seq_data%>%
  select(OTU,kingdom:species, guild)%>%
  distinct()

gs_dat<-read.csv('Processed_data/gs_dataset_mycorrhizae_Hiyang_3.3.25.csv')

gs_dat_genus<-gs_dat%>%filter(!is.na(genus))%>%group_by(genus)%>%summarise(mean_gs=mean(GS))

tax_gs_bag<-tax%>%
left_join(gs_dat_genus)%>%mutate(mean_gs = ifelse(is.na(mean_gs), "unknown", mean_gs))


#Soil data
tax_myco_soil<-read.csv('Processed_data/Soil_ITS_myco_tax.csv')

tax_gs_soil<-tax_myco_soil%>%rename(genus=Genus)%>%
  left_join(gs_dat_genus)%>%mutate(mean_gs = ifelse(is.na(mean_gs), "unknown", mean_gs))




write.csv(tax_gs_bag, 'Processed_data/taxa_w_genome_size_bag_data.csv',row.names = FALSE)
write.csv(tax_gs_soil, 'Processed_data/taxa_w_genome_size_soil_data.csv',row.names = FALSE)

          