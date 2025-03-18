#LOAD LIBRARARIES
library(tidyverse)
library(vegan)
library(readr)
library(readxl)



#Clean dat at location level
clean_dat <- read.csv('Processed_data/ITS_Site_clean.csv')


Guild_dat<-clean_dat%>%  
  filter(confidenceRanking %in% c("Highly Probable","Probable"))%>%
  mutate(guild2 = case_when(trophicMode == 'Saprotroph' ~ 'Saprotroph',
                            str_detect(guild, "Arbuscular Mycorrhizal") ~ 'Arbuscular Mycorrhizal',
                            str_detect(guild, "Ectomycorrhizal") ~ 'Ectomycorrhizal',
                            str_detect(guild, "Ericoid Mycorrhizal") ~ 'Ericoid Mycorrhizal',
                            str_detect(guild, "Orchid Mycorrhizal") ~ 'Orchid Mycorrhizal',
                            str_detect(guild, "Endophyte") ~ 'Endophyte',
                            trophicMode == 'Pathotroph' ~ 'Pathogen',
                            str_detect(guild, "Saprotroph") ~ 'Saprotroph',
                             trophicMode == 'Pathotroph-Saprotroph'~'Pathogen-Saprotroph',
                            guild=='Epiphyte'~'Other',
                            str_detect(guild, "Pathogen") ~ 'Pathogen',
                            str_detect(guild, "Lichen") ~ 'Other',
                            TRUE ~ 'unassigned'))


temp<-Guild_dat%>%filter(guild2=='unassigned')
rm(clean_dat)
#funguild output
Myco<-read.csv("Funguild/Myco_guilds.csv")#IDing all OTUs that are mycorrhizal 

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



dat_myco_hp<-myco_dat %>%  filter(confidenceRanking %in% c("Highly Probable"))

# Summarizing fungal OTUs
otu_summary_hp <- myco_dat %>%  filter(confidenceRanking %in% c("Highly Probable"))%>%
  group_by(OTU) %>%
  summarise(
    total_reads = sum(count),
    num_guilds = sum(str_count(guild, "-")) + 1,  # Count dashes to infer multiple guilds
    guilds = paste(unique(guild), collapse = "; "),  # List all unique guilds
    genera = paste(unique(genus), collapse = "; ")   # List all unique genera
  ) %>%
  mutate(
    category = ifelse(num_guilds > 1, "Multiple Guilds", "Single Guild")
  ) %>%
  group_by(category) %>%
  summarise(
    num_OTUs = n(),
    perc_total_reads = (sum(total_reads) / sum(dat_myco_hp$count)) * 100
  ) %>% rename_with(~ paste0(., "_HP"), -category)

print(otu_summary_hp)




otu_summary_p <- myco_dat %>%
  group_by(OTU,confidenceRanking) %>%
  summarise(
    total_reads = sum(count),
    num_guilds = sum(str_count(guild, "-")) + 1,  # Count dashes to infer multiple guilds
    guilds = paste(unique(guild), collapse = "; "),  # List all unique guilds
    genera = paste(unique(genus), collapse = "; "),   # List all unique genera
    family = paste(unique(family), collapse = "; ")   # List all unique family
    
    ) %>%
  mutate(
    category = ifelse(num_guilds > 1, "Multiple Guilds", "Single Guild")
  ) 

probable_summary<-otu_summary_p%>%
  group_by(category) %>%
  summarise(
    num_OTUs = n(),
    perc_total_reads = (sum(total_reads) / sum(myco_dat$count)) * 100
  ) %>% rename_with(~ paste0(., "_P"), -category)

probable_summary

multi_G<-otu_summary_p%>%filter(category=="Multiple Guilds")%>%
  group_by(confidenceRanking)%>%
  summarise(
    num_OTUs = n(),
    perc_total_reads = (sum(total_reads) / sum(myco_dat$count)) * 100
  )%>%  pivot_wider(names_from = confidenceRanking, values_from = c(num_OTUs, perc_total_reads))%>%
  mutate(category='Multiple Guilds')
multi_G

Multi_OTUs<-otu_summary_p%>%filter(category=="Multiple Guilds")%>%
  mutate( perc_total_reads = (sum(total_reads) / sum(myco_dat$count)) * 100)%>%
  arrange(desc(perc_total_reads))%>%select(-category)
Multi_OTUs<-Multi_OTUs[1:10,]


otu_summary_combined <- full_join(otu_summary_hp, multi_G, by = "category") 


write_xlsx(Multi_OTUs,'Bag_data_seq/Top_10_Multi.xlsx')


#############All


Guild_dat

Guild_dat_hp<-Guild_dat %>%  filter(confidenceRanking %in% c("Highly Probable"))

# Summarizing fungal OTUs
otu_summary_hp_all <- Guild_dat %>%  filter(confidenceRanking %in% c("Highly Probable"))%>%
  group_by(OTU) %>%
  summarise(
    total_reads = sum(count),
    num_guilds = sum(str_count(guild, "-")) + 1,  # Count dashes to infer multiple guilds
    guilds = paste(unique(guild), collapse = "; "),  # List all unique guilds
    genera = paste(unique(genus), collapse = "; ")   # List all unique genera
  ) %>%
  mutate(
    category = ifelse(num_guilds > 1, "Multiple Guilds", "Single Guild")
  ) %>%
  group_by(category) %>%
  summarise(
    num_OTUs = n(),
    perc_total_reads = (sum(total_reads) / sum(Guild_dat_hp$count)) * 100
  ) %>% rename_with(~ paste0(., "_HP"), -category)

print(otu_summary_hp_all)





otu_summary_p_all <- Guild_dat %>%
  group_by(OTU,confidenceRanking) %>%
  summarise(
    total_reads = sum(count),
    num_guilds = sum(str_count(guild, "-")) + 1,  # Count dashes to infer multiple guilds
    guilds = paste(unique(guild), collapse = "; "),  # List all unique guilds
    genera = paste(unique(genus), collapse = "; "),   # List all unique genera
    family = paste(unique(family), collapse = "; ")   # List all unique family
    
  ) %>%
  mutate(
    category = ifelse(num_guilds > 1, "Multiple Guilds", "Single Guild")
  ) 

probable_summary_all<-otu_summary_p_all%>%
  group_by(category) %>%
  summarise(
    num_OTUs = n(),
    perc_total_reads = (sum(total_reads) / sum(Guild_dat$count)) * 100
  ) %>% rename_with(~ paste0(., "_P"), -category)

probable_summary_all

multi_G_all<-otu_summary_p_all%>%filter(category=="Multiple Guilds")%>%
  group_by(confidenceRanking)%>%
  summarise(
    num_OTUs = n(),
    perc_total_reads = (sum(total_reads) / sum(myco_dat$count)) * 100
  )%>%  pivot_wider(names_from = confidenceRanking, values_from = c(num_OTUs, perc_total_reads))%>%
  mutate(category='Multiple Guilds')
multi_G_all

Multi_OTUs_all<-otu_summary_p_all%>%filter(category=="Multiple Guilds")%>%
  mutate( perc_total_reads = (sum(total_reads) / sum(myco_dat$count)) * 100)%>%
  arrange(desc(perc_total_reads))%>%select(-category)
Multi_OTUs_all<-Multi_OTUs_all[1:10,]


#otu_summary_combined <- full_join(otu_summary_hp, multi_G, by = "category") 


write_xlsx(Multi_OTUs_all,'Bag_data_seq/Top_10_Multi_all.xlsx')



#calc percentage of two OTUs

# Calculate total reads per barcode
total_reads_per_barcode <- Guild_dat %>%
  group_by(barcode) %>%
  summarise(total_reads = sum(count, na.rm = TRUE), .groups = "drop")

# Filter for the two specific OTUs and calculate their read percentages
otu_percentages <- Guild_dat %>%
  filter(OTU %in% c("ITSall_OTUk_1642", "ITSall_OTUj_1351")) %>%
  group_by(barcode, OTU) %>%
  summarise(otu_reads = sum(count, na.rm = TRUE), .groups = "drop") %>%
  left_join(total_reads_per_barcode, by = "barcode") %>%
  mutate(perc_total_reads = (otu_reads / total_reads) * 100)%>%
  left_join(Guild_dat)

p <- ggplot(otu_percentages, aes(x = barcode, y = perc_total_reads, fill = family, text = paste("Site:", Site, "<br>Transect:", Transect))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Barcode", y = "Percent of Total Reads", title = "% of Reads by sample for Selected OTUs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Convert to interactive plot with Plotly
plotly::ggplotly(p, tooltip = "text")




























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

rm(Myco,Fun_Traits)

#All meta data from 12 sites with bags collected
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')%>%
  filter(!Tube_ID %in% c(84))%>% mutate(Tube_ID=as.factor(Tube_ID))
  #Tube 84 is linked to two locations 11-1-16 and 11-2-16 due to extraction error (good for Site analysis bad for transect)
  # setdiff(names(Bag_data), names(Bag_Site_ID))
  # setdiff(names(Bag_Site_ID), names(Bag_data))


#Sample site 56 and sample site 49 both not enough biomass during harvest

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


