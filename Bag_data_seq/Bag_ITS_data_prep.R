#LOAD LIBRARARIES
library(tidyverse)
library(vegan)
library(readxl)
library(readr)
library(ggplot2)
library(stringr)

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
  #but for transect analyses they should be removed when looking at location based correlations
  mutate(Tube_ID = case_when(
    Location == 3 & Transect == 2 & Site == 11 ~ "95",
    Location == 33 & Transect == 2 & Site == 11 ~ "96",
    TRUE ~ Tube_ID  # Keep the original value if no condition is met
  ))%>%
  filter(!Tube_ID %in% c('84'))%>%#Tube 84 is linked to two locations 11-1-16 and 11-2-16 due to extraction error (good for Site analysis bad for transect)
  filter(!Tube_ID %in% c(8,11))#these tubes are both from site 29, but I mixed up the transects and it creates duplication errors when I keep samples in
  rm(CNP_dat,dat)

Guild_dat<-clean_dat%>%  filter(confidenceRanking %in% c('Probable',"Highly Probable"))%>%
  mutate(guild2 = case_when(trophicMode == 'Saprotroph' ~ 'Saprotroph',
                            str_detect(guild, "Arbuscular Mycorrhizal") ~ 'Arbuscular Mycorrhizal',
                            str_detect(guild, "Ectomycorrhizal") ~ 'Ectomycorrhizal',
                            str_detect(guild, "Orchid Mycorrhizal") ~ 'Orchid Mycorrhizal',
                            str_detect(guild, "Ericoid Mycorrhizal") ~ 'Ericoid Mycorrhizal',
                            str_detect(guild, "Endophyte") ~ 'Endophyte',
                            str_detect(guild, "Epiphyte") ~ 'Epiphyte',
                            str_detect(family, "Mortierellaceae") ~ 'Saprotroph',
                            str_detect(guild, "Parasite") ~ 'Parasite',
                            str_detect(guild, "Pathogen") ~ 'Pathogen',
                            TRUE ~ 'unassigned'))

temp<-Guild_dat%>%filter(guild2=='unassigned')
rm(clean_dat)
#funguild output
Myco<-read.csv("Funguild/Myco_guilds.csv")#IDing all OTUs that are mycorrhizal 

#explo types Polme et al 2021
Fun_Traits<-read_excel("Processed_data/Polme_etal_2021_FungalTraits.xlsx")

#select only mycorrhizal taxa
myco_dat<-Guild_dat%>%
  mutate (guild = case_when(genus %in% Myco$Genus ~ 'mycorrhizal',
                            TRUE ~ guild ))%>%
  filter(guild=='mycorrhizal')%>%
  left_join(Fun_Traits%>%rename(genus=GENUS)%>%
              select(genus,Ectomycorrhiza_exploration_type_template,Ectomycorrhiza_lineage_template))%>%
  rename(exploration_type=Ectomycorrhiza_exploration_type_template,
         Ecm_lineage=Ectomycorrhiza_lineage_template)%>%
  mutate(across(c(exploration_type,Ecm_lineage), ~replace_na(.x, "unknown")))# Replace NA with "Unknown" 

#select taxa
all_tax<-Guild_dat%>%
  select(OTU,kingdom:species,confidenceRanking, guild,guild2)%>%
  distinct()

myco_tax<-myco_dat%>%
  select(OTU,kingdom:species, guild,guild2)%>%
  distinct()

rm(Myco,Fun_Traits)

#All meta data from 12 sites with bags collected
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')[,-1]

Bag_Site_ID<-read.csv('Processed_data/All_Bag_Site_Info.csv')%>%
  filter(!Tube_ID %in% c(8,11))#these tubes are both from site 29, but I mixed up the transects and 
  #it creates duplication errors when I keep samples in
  #I could include them in Site level analyses.
  # setdiff(names(Bag_data), names(Bag_Site_ID))
  # setdiff(names(Bag_Site_ID), names(Bag_data))

Bag_data<-Bag_data%>%left_join(Bag_Site_ID, relationship = "many-to-many")%>%#doing this to get Tube_ID from Bag_Site_ID
  filter(!is.na(Tube_ID))%>%mutate(Tube_ID=as.factor(Tube_ID))#remove cols with na values left over from removing tubes 8 and 11

rm(Bag_Site_ID)
#t<-Bag_data%>%select(Site,Transect,Location,Tube_ID)%>%group_by(Site,Transect)%>% mutate(n())
#missing 2 samples from site 29 crossing transects
#Sample site 56 and sample site 49

#All taxa
Bag_Guild<-Bag_data%>%
  dplyr::select(Tube_ID,Site,Transect,Fire.Severity,Fire.Interval, Location)%>%
  right_join(Guild_dat%>% mutate(Tube_ID=as.factor(Tube_ID)))

#Myco taxa
Bag_myco<-Bag_data%>%
  dplyr::select(Tube_ID,Site,Transect,Fire.Severity,Fire.Interval, Location)%>%
  right_join(myco_dat%>% mutate(Tube_ID=as.factor(Tube_ID)))

#Guild Comp per Location
Bag_Guild %>%
  group_by(site_transect_locat = paste("Site", Site, "T", Transect, 'L',Location)) %>%
  mutate(
    Site = as.factor(Site),
    guild2 = factor(guild2, levels = c(setdiff(unique(guild2), "other_unknown"), "other_unknown")),
    # Calculate relative abundance within each site/transect/loc
    rel_abundance = count/sum(count) * 100,
    # Create a factor with custom ordering for site/transect combinations
    site_transect = factor(site_transect_locat,
                           levels = unique(site_transect_locat)[order(Site, Transect,Location)]))%>%
  ggplot(aes(x = site_transect_locat, y = rel_abundance, fill = guild2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Fire.Interval, scales = "free_x") +
  labs(
    x = "Site/Transect",
    y = "Relative Abundance (%)",
    fill = "Guild" ) 


#create sample rarefecation curve
temp <- Bag_Guild %>%
  select(Tube_ID, OTU, count) %>%
  pivot_wider(names_from = OTU, values_from = count, values_fill = 0)

mat<-temp%>%
  remove_rownames()%>%
  column_to_rownames("Tube_ID")

# assess variation in sampling effort, plotting sample effort curves
#mat <- mat[rowSums(mat) >= 1000, ]  # Keep only samples with at least 1000 reads
curve<- rarecurve(mat, step=1000, tidy=TRUE)

Bag_data%>%mutate(Tube_ID=as.factor(Tube_ID)) %>%
  left_join( curve%>%rename(Tube_ID=Site,Samples=Sample),by = join_by('Tube_ID'))%>% 
  ggplot(aes(x=Samples, y=Species, colour=as.factor(Sample), group=Tube_ID)) + 
  geom_line()




#format required for some vegan functions
wide_myco<-myco_dat  %>%
  pivot_wider(
    names_from = OTU,         # Column containing OTU names that will become new columns
    values_from = count,      # Column containing values that will fill the new columns
    id_cols = Tube_ID,  # Column(s) to keep as identifier
    values_fill = 0
  )%>% mutate(Tube_ID=as.factor(Tube_ID))

Bag_Seq_wide<-left_join(Bag_data%>%
                          select(Tube_ID,Fire.Interval,Fire.Severity,Site, Transect, Location,#Site Data
                                 Ortho_P_mg_kg,Nitrate_mg_kg, Ammonia_mg_kg,#Nutrients from resin data
                                 log10_biomass_day,perc_myco_host_freq, #biomass data
                                 Longitude,Latitude), #meta data
                        wide_myco)%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))

#check to make sure there are no NAs
# NA_rows <- Bag_Seq_wide %>%
#   mutate(across(everything(), ~ ifelse(is.na(.), ., NA))) %>%
#   filter(if_any(everything(), ~ !is.na(.)))

write.csv(Bag_Seq_wide, 'Processed_data/Bag_Seq_wide.csv', row.names = FALSE)
write.csv(Bag_data, 'Processed_data/Updated_Bag_data.csv', row.names = FALSE)
write.csv(myco_tax, 'Processed_data/Bag_dat_myco_tax.csv', row.names = FALSE)


#myco
myco_reads<-Bag_myco%>%
  summarise(sum(count))%>%
  pull()


dat_myco_RA<-Bag_myco%>%
  group_by(Site,Transect,Location)%>%
  summarise(reads_samp=sum(count))%>%
  left_join(Bag_myco)%>%
  left_join(Bag_myco%>%
              group_by(Fire.Severity)%>%
              summarise(severity_reads=sum(count)))%>%
  left_join(Bag_myco%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads=sum(count)))%>%
  mutate(RA_samp= count/reads_samp,
         RA_total_reads= count/myco_reads,
         RA_total_interval= count/interval_reads,
         RA_total_severity= count/severity_reads)%>%
  left_join(Bag_data%>%select(Site,Transect,Location,perc_myco_host_freq,
                              Nitrate_mg_kg,Ammonia_mg_kg,Ortho_P_mg_kg))
#all seq data from bags
total_reads<-Bag_Guild%>%
  summarise(sum(count))%>%
  pull()


dat_all_RA<-Bag_Guild%>%
  group_by(Site,Transect,Location)%>%
  summarise(reads_samp=sum(count))%>%
  left_join(Bag_Guild)%>%
  left_join(Bag_Guild%>%
              group_by(Fire.Severity)%>%
              summarise(severity_reads=sum(count)))%>%
  left_join(Bag_Guild%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads=sum(count)))%>%
  mutate(RA_samp= count/reads_samp,
         RA_total_reads= count/total_reads,
         RA_total_interval= count/interval_reads,
         RA_total_severity= count/severity_reads)

