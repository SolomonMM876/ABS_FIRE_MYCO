library(tidyverse)

Summary<-read_tsv('Processed_data/CG/ITS_data_quality_summary.tsv')
TEmp<-read_tsv('Processed_data/CG/ITS_output_clean.tsv')


description<-Summary%>%
  #this removes the samples that were not properly sequenced (Sites 49 Tran 2- Site 63)
  filter(!str_detect(barcode, "^CG124\\-ITS\\-CG(9[6-9]|1[0-1][0-9]|12[0-4])$"))


Decomp_Site<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Updated_Soil_data.csv')

Decomp_Site%>%select(barcode)%>%distinct()%>%arrange(desc(barcode))->t
 #49 samples to start

Guild_dat <- Decomp_Site %>%
  filter(confidenceRanking == "Highly Probable" | OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073")) %>%  # Keep "Highly Probable" & specified OTUs
  filter(!str_detect(guild, "-") | OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073")) %>%  # Exclude multiple guilds except for specified OTUs
  mutate(
    guild2 = case_when(
      str_detect(guild, "Arbuscular Mycorrhizal") ~ "Arbuscular Mycorrhizal",
      str_detect(guild, "Ectomycorrhizal") ~ "Ectomycorrhizal",
      str_detect(guild, "Ericoid Mycorrhizal") ~ "Ericoid Mycorrhizal",
      str_detect(guild, "Orchid Mycorrhizal") ~ "Orchid Mycorrhizal",
      str_detect(guild, "Saprotroph") ~ 'Saprotroph',
      trophicMode == 'Pathotroph' ~ 'Pathogen',
      str_detect(guild, "Endophyte") ~ 'Endophyte',
      guild=='Epiphyte'~'Other',
      str_detect(guild, "Lichen") ~ 'Other'
    ),
    genus = case_when(  # Separate case_when() for genus assignment
      OTU == 'ITSall_OTUa_10308' ~ 'Tomentella',
      OTU == 'ITSall_OTUk_1642' ~ 'Coltricia',
      OTU == 'ITSall_OTUm_3073' ~ 'Ionosporus',
      TRUE ~ genus  # Keep the existing genus if OTU is not listed
    ),
    note = if_else(OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073"), 
                   "Included due to manual BLAST confirmation", 
                   NA_character_))  # Add note for manually included OTUs
Guild_dat%>%select(barcode)%>%distinct()%>%arrange(desc(barcode))->t


#explo types Polme et al 2021
Fun_Traits<-read_excel("Processed_data/Polme_etal_2021_FungalTraits.xlsx")

#select only mycorrhizal taxa
myco_dat<-Guild_dat%>%
  filter(str_detect(guild, "mycorrhizal") | str_detect(guild, "Mycorrhizal")) %>%
  left_join(Fun_Traits%>%rename(genus=GENUS)%>%
              select(genus,Ectomycorrhiza_exploration_type_template,Ectomycorrhiza_lineage_template))%>%
  rename(exploration_type=Ectomycorrhiza_exploration_type_template,
         Ecm_lineage=Ectomycorrhiza_lineage_template)%>%
  mutate(exploration_type = na_if(exploration_type, "unknown"))  # Convert "Unknown" to NA


myco_dat<-myco_dat%>% #add the total_myco reads per sample
  left_join(myco_dat%>%group_by(barcode)%>%summarise(myco_reads=sum(count)))%>%
  mutate(barcode=as.factor(barcode))

Decomp_Site %>% summarise(sum(count,na.rm=T))
myco_dat %>% summarise(sum(count, na.rm=T))

myco_dat%>%select(barcode)%>%distinct()%>%arrange(desc(as.numeric(barcode)))->t

#Save myco data

write.csv(myco_dat, 'Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_myco_dat.csv', row.names = FALSE)

#make OTU table
otu_table<-myco_dat%>%select(barcode, OTU, count) %>% 
  pivot_wider(names_from = OTU, values_from = count, values_fill = 0)%>%
  column_to_rownames('barcode')


write.csv(otu_table,'Soil_ITS_scripts/Jeff_proc/Processed_data/otu_table_soil.csv')



#select taxa
all_tax_soil<-Guild_dat%>%
  select(OTU,kingdom:species,confidenceRanking, guild,guild2)%>%
  distinct()

myco_tax<-myco_dat%>%
  select(OTU,kingdom:species, guild,guild2)%>%
  distinct()

#Redo with just selected sites

#All 
Soil_Guild<-Guild_dat%>% 
  dplyr::select(barcode:Transect,Fire.Severity,Fire.Interval)

#Myco
Soil_myco<-myco_dat%>%
  dplyr::select(barcode:Transect,Fire.Severity,Fire.Interval,exploration_type)


#format required for some vegan functions
wide_myco<-myco_dat  %>%
  pivot_wider(
    names_from = OTU,         # Column containing OTU names that will become new columns
    values_from = count,      # Column containing values that will fill the new columns
    id_cols = barcode,  # Column(s) to keep as identifier
    values_fill = 0
  )%>% mutate(barcode=as.factor(barcode))

Soil_Seq_wide<-right_join(Guild_dat %>% select(barcode,sample:Dead.Tree.Canopy.Cover_perc) %>% distinct()
                          ,wide_myco)%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))



#Save OUTputs


write.csv(Soil_Seq_wide, 'Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_Seq_wide.csv', row.names = FALSE)
write.csv(Decomp_Site, 'Soil_ITS_scripts/Jeff_proc/Processed_data/Updated_Soil_data.csv', row.names = FALSE)
write.csv(myco_tax, 'Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_dat_myco_tax.csv', row.names = FALSE)


#myco


dat_myco_RA_soil<-Soil_myco%>%
  group_by(Site,Transect)%>%
  summarise(reads_samp=sum(count))%>%
  left_join(Soil_myco)%>%
  left_join(Soil_myco%>%
              group_by(Fire.Severity)%>%
              summarise(severity_reads=sum(count)))%>%
  left_join(Soil_myco%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads=sum(count)))%>%
  ungroup() %>% 
  mutate(myco_reads=sum(count),
         RA_samp= count/reads_samp,
         RA_total_reads= count/myco_reads,
         RA_total_interval= count/interval_reads,
         RA_total_severity= count/severity_reads)

dat_explo_soil<-dat_myco_RA_soil%>%
  filter(!is.na(exploration_type))%>%
  group_by(Site,Transect,exploration_type)%>%
  summarise(explo_count=sum(count))%>%
  left_join(dat_myco_RA_soil%>%select(Site,Transect,barcode,Fire.Severity,Fire.Interval,reads_samp)%>%distinct())%>%
  left_join(dat_myco_RA_soil%>%filter(!is.na(exploration_type))%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads_explo=sum(count)))%>%
  left_join(dat_myco_RA_soil%>%filter(!is.na(exploration_type))%>%
              group_by(Fire.Severity)%>%
              summarise(Severity_reads_explo=sum(count)))%>%
  mutate(RA_explo_Interval=explo_count/interval_reads_explo,
         log_RA_explo_Interval=log10(RA_explo_Interval+(1.490598e-05)/2),
         RA_explo_Severity=explo_count/Severity_reads_explo,
         log_RA_explo_Severity=log10(RA_explo_Severity+(1.451303e-05)/2))


#all seq data from bags

#this is a more broad guild alignment from include script
Guild_dat_hp<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_data_all_guild.csv')

Soil_Guild<-Guild_dat_hp%>% 
  dplyr::select(barcode:Transect,Fire.Severity,Fire.Interval)



dat_all_RA_soil<-Soil_Guild%>%
  group_by(Site,Transect,barcode)%>%
  summarise(reads_samp=sum(count))%>%
  left_join(Soil_Guild)%>%
  ungroup()%>%
  left_join(Soil_Guild%>%
              group_by(Fire.Severity)%>%
              summarise(severity_reads=sum(count)))%>%
  ungroup()%>%
  left_join(Soil_Guild%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads=sum(count)))%>%
  ungroup()%>%
  mutate(total_reads=sum(count),
         RA_samp= count/reads_samp,
         RA_total_reads= count/total_reads,
         RA_total_interval= count/interval_reads,
         RA_total_severity= count/severity_reads)














