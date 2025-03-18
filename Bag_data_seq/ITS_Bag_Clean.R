library(tidyverse)

dat <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_clean.tsv')

#filter out data for CNP measurements and community
CNP_dat<-dat%>%
  filter(is.na(Site))

write.csv(CNP_dat, file='Processed_data/CNP_seq_dat.csv',row.names=FALSE)

#Clean dat at location level
clean_dat <- dat %>%
  #remove transect pooled samples for CN
  anti_join(CNP_dat, by = colnames(CNP_dat)) %>%
  #Samples 95 and 96 were pooled. I have separated them,
  # for transect analyses they should be removed when looking at location based correlations
  mutate(Tube_ID = case_when(
    Location == 3 & Transect == 2 & Site == 11 ~ "95",
    Location == 33 & Transect == 2 & Site == 11 ~ "96",
    TRUE ~ Tube_ID),  # Keep the original value if no condition is met
    barcode = case_when(
      Location == 3 & Transect == 2 & Site == 11 ~ "SMM95",
      Location == 33 & Transect == 2 & Site == 11 ~ "SMM96",
      TRUE ~ barcode)
  )%>%
  filter(!Tube_ID %in% c('84'))#Tube 84 is linked to two locations 11-1-16 and 11-2-16 due to extraction error (good for Site analysis bad for transect)

Tube_811<-clean_dat%>%
  filter(Tube_ID %in% c(8,11))%>%#these tubes are both from the same location, but I processed two samples here so I am combining them and added it back as one sample
  group_by(OTU, kingdom, phylum, class, order, family, genus, species, SH_species, SH_gb_acc, SH_sh_acc, guild, confidenceRanking, 
           growthForm,trait, Site, Transect, Location) %>%
  summarise(
    count = round(as.numeric(mean(count, na.rm = TRUE)),digits = 0),
    resampled_count = round(as.numeric(mean(resampled_count, na.rm = TRUE)),digits = 0),
  )%>%
  mutate(Tube_ID= '8',
         barcode= 'SMM8')

clean_dat<-clean_dat%>%
  #remove tubes before adding summarized rows
  filter(!Tube_ID %in% c(8,11))

clean_dat<-bind_rows(clean_dat,Tube_811)
  

rm(CNP_dat,dat,Tube_811)

clean_dat%>%select(Tube_ID)%>%distinct()%>%arrange(desc(as.numeric(Tube_ID)))->t


#I have removed Tube ID 11 from this df as well as Tube ID 84
write.csv(clean_dat, 'Processed_data/ITS_Site_clean.csv',row.names = FALSE)
