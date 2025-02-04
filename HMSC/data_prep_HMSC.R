library(readr)
library(tidyr)
library(dplyr)


#combined bag data
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')%>%
  mutate(Location_ID = paste0("S", Site, "T", Transect, "L", Location))
#bag_data with tubeIDS
Bag_Site_ID<-read.csv('Processed_data/All_Bag_Site_Info.csv')%>%
  filter(!Tube_ID %in% c(8,11))#these tubes are both from site 29, but I mixed up the transects and 
#it creates duplication errors when I keep samples in
#I could include them in Site level analyses.

Bag_data<-Bag_data%>%left_join(Bag_Site_ID, relationship = "many-to-many")
rm(Bag_Site_ID)

#read in community data
Community<-read.csv('Processed_data/cleaned_seq_dat.csv')%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))%>%
  filter(guild %in%c("Ectomycorrhizal","Arbuscular Mycorrhizal","Orchid Mycorrhizal","Ericoid Mycorrhizal"))

dat_summary <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_summarised.tsv')%>%
  left_join(Community%>%select(Tube_ID,barcode)%>%distinct())%>%filter(!is.na(Tube_ID))

# Data prep step
data <- Community %>%
  filter(!Tube_ID%in% c(8,11))%>% #these samples are excluded because they are same location, but different transects
  group_by(Tube_ID, OTU) %>%
  #transform to wide format
  summarise(sequence_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(values_from = sequence_count, names_from = OTU,  values_fill = 0) %>%
  left_join(Bag_data%>% 
              distinct(Tube_ID,Site,Transect,Location,Fire.Interval, Fire.Severity,
                       log10_biomass_day,C_N,C_P,N_P,
                       Ortho_P_mg_kg,Nitrate_mg_kg,Ammonia_mg_kg,pH,
                       perc_myco_host_freq,
                       Latitude,Longitude),
            by = "Tube_ID")%>%
  left_join(dat_summary%>% select(Tube_ID,n_reads))%>%
  rename(readcount=n_reads)%>%
  mutate(Tube_ID=as.factor(Tube_ID),
         Site=as.factor(Site),
         Transect=as.factor(Transect),
         Location=as.factor(Location))%>%
  #left_join(Bag_data%>%select(Location_ID, Tube_ID), by= 'Tube_ID')%>%
  select(Tube_ID,Site,Transect,Location, Fire.Interval, Fire.Severity,readcount,
         log10_biomass_day,C_N,N_P,N_P,C_P,
         perc_myco_host_freq, 
         Latitude,Longitude,
         Ortho_P_mg_kg,Nitrate_mg_kg,Ammonia_mg_kg,pH,
         everything())#


rm(dat_summary,Bag_data)



TrData<-Community %>%
  filter(!Tube_ID%in% c(8,11))%>% #these samples are excluded because they are same location, but different transects
  select( kingdom:species,SH_species,OTU) %>% distinct()%>%
  mutate(across(everything(), ~ ifelse(is.na(.), "Unknown", .))) # Replace NA with "Unknown" 
         # I have to do this in order to make the phylo tree later, it is not a perfect solution, but it is the best I can think of
         
library(readxl)
Fun_Traits<-read_excel("Processed_data/Polme_etal_2021_FungalTraits.xlsx")
# 
 TrData<-TrData%>%
   left_join(Fun_Traits, by = c('genus'='GENUS'))%>%
   select(kingdom:OTU,Ectomycorrhiza_exploration_type_template,Ectomycorrhiza_lineage_template)%>%
   rename(exploration_type=Ectomycorrhiza_exploration_type_template,
          Ecm_lineage=Ectomycorrhiza_lineage_template)%>%
   mutate(across(everything(), ~replace_na(.x, "unknown")))

rm(Community,Fun_Traits)


write.csv(data, file='HMSC/data/Bag_data.csv',row.names=FALSE)
write.csv(TrData, file='HMSC/data/Trait_Phylo_data.csv',row.names=FALSE)

