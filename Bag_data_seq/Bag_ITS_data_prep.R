#LOAD LIBRARARIES
library(tidyverse)
library(vegan)
library(readr)
library(readxl)



#Clean dat at location level
clean_dat <- read.csv('Processed_data/ITS_Site_clean.csv')

clean_dat%>%select(Tube_ID)%>%distinct()%>%arrange(desc(as.numeric(Tube_ID)))->t
#90 samples to start


Guild_dat <- clean_dat %>%
  filter(confidenceRanking == "Highly Probable" | OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073")) %>%  # Keep "Highly Probable" & specified OTUs
  filter(!str_detect(guild, "-") | OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073")) %>%  # Exclude multiple guilds except for specified OTUs
  mutate(guild2 = case_when(
    str_detect(guild, "Arbuscular Mycorrhizal") ~ "Arbuscular Mycorrhizal",
    str_detect(guild, "Ectomycorrhizal") ~ "Ectomycorrhizal",
    str_detect(guild, "Ericoid Mycorrhizal") ~ "Ericoid Mycorrhizal",
    str_detect(guild, "Orchid Mycorrhizal") ~ "Orchid Mycorrhizal",
    TRUE ~ "other"
  ),
  note = if_else(OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073"), 
                 "Included due to manual BLAST confirmation", 
                 NA_character_))  # Add note for manually included OTUs

Guild_dat%>%select(Tube_ID)%>%distinct()%>%arrange(desc(as.numeric(Tube_ID)))->t
#lost one sample in filtering for taxa with only highly probable and single guild

#temp<-Guild_dat%>%filter(guild2=='Ericoid Mycorrhizal')
rm(clean_dat)


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
  mutate(Tube_ID=as.factor(Tube_ID))

myco_dat%>%select(Tube_ID)%>%distinct()%>%arrange(desc(as.numeric(Tube_ID)))->t
#lost 2 samples filtering for mycorrhizal OTUs

write.csv(myco_dat, 'Processed_data/Bag_Seq_myco_dat.csv', row.names = FALSE)




otu_table<-myco_dat%>%select(Tube_ID, OTU, count) %>% 
  pivot_wider(names_from = OTU, values_from = count, values_fill = 0)%>%
  column_to_rownames('Tube_ID')#Tube_ID 77 and 89 have no mycorrhizal fungi in them 
  

write.csv(otu_table,'Processed_data/otu_table_bag.csv')

#select taxa
all_tax<-Guild_dat%>%
  select(OTU,kingdom:species,confidenceRanking, guild,guild2)%>%
  distinct()

myco_tax<-myco_dat%>%
  select(OTU,kingdom:species, guild,guild2)%>%
  distinct()

rm(Fun_Traits)


#All meta data from 12 sites with bags collected
#from df joins script
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')%>%
  filter(!Tube_ID %in% c(84))%>% mutate(Tube_ID=as.factor(Tube_ID))
  #Tube 84 is linked to two locations 11-1-16 and 11-2-16 due to extraction error (good for Site analysis bad for transect)
  # setdiff(names(Bag_data), names(Bag_Site_ID))
  # setdiff(names(Bag_Site_ID), names(Bag_data))


#1 Sample site 56 and 1 sample site 49 both not enough biomass during harvest

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
                                 Ortho_P_mg_kg,Nitrate_mg_kg, Ammonia_mg_kg,pH,#Nutrients from resin data
                                 log10_biomass_day,perc_myco_host_freq, #biomass data
                                 C_N,C_P,N_P,Carb_Hyph,Nitrog_Hyph,Phos_Hyph,
                                 Longitude,Latitude)%>%
                          left_join(myco_dat%>%select(Tube_ID,myco_reads)%>%distinct()), #meta data
                        wide_myco)%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))

#check to make sure there are no NAs
# NA_rows <- Bag_Seq_wide %>%
#   mutate(across(everything(), ~ ifelse(is.na(.), ., NA))) %>%
#   filter(if_any(everything(), ~ !is.na(.)))

write.csv(Bag_Seq_wide, 'Processed_data/Bag_Seq_wide.csv', row.names = FALSE)#bag data with selected meta and Myco communities
write.csv(Bag_data, 'Processed_data/Updated_Bag_data.csv', row.names = FALSE)#Just Site info, veg nute meta etc
write.csv(myco_tax, 'Processed_data/Bag_dat_myco_tax.csv', row.names = FALSE)#mycorrhizal taxa included


#myco


dat_myco_RA_bag<-Bag_myco%>%
  group_by(Site,Transect,Location)%>%
  summarise(reads_samp=sum(count))%>%
  left_join(Bag_myco)%>%
  left_join(Bag_myco%>%
              group_by(Fire.Severity)%>%
              summarise(severity_reads=sum(count)))%>%
  left_join(Bag_myco%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads=sum(count)))%>%
  mutate(myco_reads=sum(count),
    RA_samp= count/reads_samp,
         RA_total_reads= count/myco_reads,
         RA_total_interval= count/interval_reads,
         RA_total_severity= count/severity_reads)%>%
  left_join(Bag_data%>%select(Site,Transect,Location,perc_myco_host_freq,
                              Nitrate_mg_kg,Ammonia_mg_kg,Ortho_P_mg_kg))

dat_explo_bag<-dat_myco_RA_bag%>%
  filter(!is.na(exploration_type))%>%
  group_by(Site,Transect,Location,exploration_type)%>%
  summarise(explo_count=sum(count))%>%
  left_join(dat_myco_RA_bag%>%select(Site,Transect,Location,Fire.Severity,Fire.Interval,reads_samp)%>%distinct())%>%
  left_join(dat_myco_RA_bag%>%filter(!is.na(exploration_type))%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads_explo=sum(count)))%>%
  left_join(dat_myco_RA_bag%>%filter(!is.na(exploration_type))%>%
              group_by(Fire.Severity)%>%
              summarise(Severity_reads_explo=sum(count)))%>%
  mutate(RA_explo_Interval=explo_count/interval_reads_explo,
         log_RA_explo_Interval=log10(RA_explo_Interval+(1.490598e-05)/2),
         RA_explo_Severity=explo_count/Severity_reads_explo,
         log_RA_explo_Severity=log10(RA_explo_Severity+(1.451303e-05)/2))


#all seq data from bags

dat_all_RA_bag<-Bag_Guild%>%
  group_by(Site,Transect,Location)%>%
  summarise(reads_samp=sum(count))%>%
  left_join(Bag_Guild)%>%
  ungroup()%>%
  left_join(Bag_Guild%>%
              group_by(Fire.Severity)%>%
              summarise(severity_reads=sum(count)))%>%
  ungroup()%>%
  left_join(Bag_Guild%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads=sum(count)))%>%
  ungroup()%>%
  mutate(total_reads=sum(count),
    RA_samp= count/reads_samp,
         RA_total_reads= count/total_reads,
         RA_total_interval= count/interval_reads,
         RA_total_severity= count/severity_reads)


