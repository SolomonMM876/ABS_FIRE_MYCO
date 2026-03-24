#LOAD LIBRARARIES
library(tidyverse)
library(vegan)
library(readr)
library(readxl)
library(writexl)


#import Blast ID df and clean data
Blast_ID<-read_excel('Raw_data/Updated_Data/ABS.MER.fielddata.Feb.2023_Site.Info_AF.xlsx')
Blast_ID$Site= sub(c('ABS00|ABS0'),'',Blast_ID$Site)
Blast_ID$sample_ID<-as.character(Blast_ID$sample_ID)
Blast_ID$transect<-as.character(Blast_ID$transect)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
names(Blast_ID)[9]<-'Veg_Class_Abv'





#Nute and Veg data
Nutrients_Transects<-read_excel('Processed_data/Nutrients_Transect_level.xlsx')
VEG_COVER_Transects <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                  sheet = "Transect.Level_Data")
VEG_COVER_Transects$Site= sub(c('ABS00|ABS0'),'',VEG_COVER_Transects$Site)
VEG_COVER_Transects$Transect= sub(c('T'),'',VEG_COVER_Transects$Transect)
Myco_plant_spp<-read.csv('Processed_data/Myco_host_abundance.csv')%>%
  mutate(Site=as.factor(Site),
         Transect=as.factor(Transect))


Blast_ID<-Blast_ID%>%
  dplyr::mutate(Regime = paste(Interval, Severity, sep = "_"))%>%
  rename(Transect=transect)%>%
  select(-`Fire 3`,-`Interval (yrs)...16`,-`FESM severity category`)%>%
  left_join(Myco_plant_spp)%>%
  left_join(Nutrients_Transects)%>%
  left_join(VEG_COVER_Transects)%>%
  #this is selecting for the 30 decomp sites which I have CN data for all of
  filter(!is.na(Carbon)) %>% 
  rename(barcode=sample_ID) %>% 
  mutate(barcode=str_replace_all(barcode,'\\.','-'))

summary(Blast_ID)#no na;s in df

Blast_ID %>% distinct(Site,Severity,Interval,Regime) %>% group_by(Severity,Interval) %>% 
  summarise(n())

Blast_ID%>% distinct(Site,Severity,Interval,Regime) %>% group_by(Regime) %>% 
  summarise(n())


#Clean dat at location level
clean_dat <- read_tsv('Processed_data/CG/ITS_output_clean.tsv')


Guild_dat<-clean_dat%>%  select(-c(sample,weight_mg,Sample_Type,qbit_DNA_conc_ng_uL)) %>% 
  filter(confidenceRanking %in% c("Highly Probable","Probable"))%>%
  mutate(guild2 = case_when(trophicMode == 'Saprotroph' ~ 'Saprotroph',
                            str_detect(guild, "Arbuscular Mycorrhizal") ~ 'Arbuscular Mycorrhizal',
                            str_detect(guild, "Ectomycorrhizal") ~ 'Ectomycorrhizal',
                            str_detect(guild, "Ericoid Mycorrhizal") ~ 'Ericoid Mycorrhizal',
                            str_detect(guild, "Orchid Mycorrhizal") ~ 'Orchid Mycorrhizal',
                            str_detect(guild, "Saprotroph") ~ 'Saprotroph',
                            trophicMode == 'Pathotroph' ~ 'Pathogen',
                            str_detect(guild, "Endophyte") ~ 'Endophyte',
                            trophicMode == 'Pathotroph-Saprotroph'~'Pathogen-Saprotroph',
                            guild=='Epiphyte'~'Other',
                            str_detect(guild, "Pathogen") ~ 'Pathogen',
                            str_detect(guild, "Lichen") ~ 'Other',
                            TRUE ~ 'unassigned'))


temp<-Guild_dat%>%filter(guild2=='unassigned')
#rm(clean_dat)


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

#One OTU is probable and is present @>5%:
#>ITSall_OTUa_5074
#aacgcatcttgcgctccttggtattccgaggagcatgcctgtttgagtgtcgtgaagagaggagtcgaacgccttgtgcaaaaggcgttcggatttggacgctgtcggaccttggtccgactcgtctcgaaatgtattagccgtcacccaaacctcatagacggacggtgtgataagttgatcgcccgagtctcgcgtgggtagggtcggcttctagtcgtcttcggacaatcagtcgttctgacctcaaatcaggtaggactacccgctgaacttaa
#UNABLE TO CONFIRM -> EXCLUDING
#see MEthods/Myco_OTUs for blast


otu_summary_combined <- full_join(otu_summary_hp, multi_G, by = "category") 


write_xlsx(Multi_OTUs,'Soil_ITS_scripts/Jeff_proc/Top_10_Multi.xlsx')


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
    perc_total_reads = (sum(total_reads) / sum(Guild_dat$count)) * 100
  )%>%  pivot_wider(names_from = confidenceRanking, values_from = c(num_OTUs, perc_total_reads))%>%
  mutate(category='Multiple Guilds')
multi_G_all

Multi_OTUs_all<-otu_summary_p_all%>%filter(category=="Multiple Guilds")%>%
  mutate( perc_total_reads = (sum(total_reads) / sum(Guild_dat$count)) * 100)%>%
  arrange(desc(perc_total_reads))%>%select(-category)
Multi_OTUs_all<-Multi_OTUs_all[1:10,]

Multi_OTUs_all

#otu_summary_combined <- full_join(otu_summary_hp, multi_G, by = "category") 


write_xlsx(Multi_OTUs_all,'Soil_ITS_scripts/Jeff_proc/top_10_Multi_all.xlsx')



#calc percentage of two OTUs

# Calculate total reads per barcode
total_reads_per_barcode <- Guild_dat %>%
  group_by(barcode) %>%
  summarise(total_reads = sum(count, na.rm = TRUE), .groups = "drop")

# Filter for the two specific OTUs and calculate their read percentages
otu_percentages <- Guild_dat %>%
  group_by(barcode, OTU) %>%
  summarise(otu_reads = sum(count, na.rm = TRUE), .groups = "drop") %>%
  left_join(total_reads_per_barcode, by = "barcode") %>%
  mutate(perc_total_reads = (otu_reads / total_reads) * 100)%>%
  left_join(Guild_dat)%>%
  left_join(Blast_ID)

p <- ggplot(otu_percentages, aes(x = barcode, y = perc_total_reads, fill = guild2, text = paste("Site:", Site, "<br>Transect:", Transect))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Barcode", y = "Percent of Total Reads", title = "% of Reads by sample for Selected OTUs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Convert to interactive plot with Plotly
plotly::ggplotly(p, tooltip = "text")






#Make rare curve

rm(Fun_Traits)

#All meta data from sites 
Decomp_Site<-Guild_dat%>%
  right_join(Blast_ID)%>%
  #this removes the samples that were not properly sequenced (Sites 49 Tran 2- Site 63)
  filter(!str_detect(barcode, "^CG124\\-ITS\\-CG(9[6-9]|1[0-1][0-9]|12[0-4])$"))


#All taxa
Soil_Guild<-Blast_ID%>% 
  dplyr::select(barcode,Site,Transect,Fire.Severity,Fire.Interval)%>%
  left_join(Guild_dat%>% mutate(barcode=as.factor(barcode)))


#Guild Comp per Location
Soil_Guild %>%
  group_by(site_transect = paste("Site", Site, "T", Transect)) %>%
  mutate(
    Site = as.factor(Site),
    #guild2 = factor(guild2, levels = c(setdiff(unique(guild2), "other_unknown"), "other_unknown")),
    # Calculate relative abundance within each site/transect/loc
    rel_abundance = count/sum(count) * 100,
    # Create a factor with custom ordering for site/transect combinations
    site_transect = factor(site_transect,
                           levels = unique(site_transect)[order(Site, Transect)]))%>%
  ggplot(aes(x = site_transect, y = rel_abundance, fill = guild2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Fire.Interval, scales = "free_x") +
  labs(
    x = "Site/Transect",
    y = "Relative Abundance (%)",
    fill = "Guild" ) 

All_OTUs<-Blast_ID%>% 
  dplyr::select(barcode,Site,Transect,Fire.Severity,Fire.Interval)%>%
  left_join(clean_dat%>% mutate(barcode=as.factor(barcode)))


#create sample rarefecation curve
temp <- All_OTUs %>%
  filter(!str_detect(barcode, "^CG124\\-ITS\\-CG(9[6-9]|1[0-1][0-9]|12[0-4])$")) %>% 
  select(barcode, OTU, count) %>%
  pivot_wider(names_from = OTU, values_from = count, values_fill = 0)

mat<-temp%>%
  remove_rownames()%>%
  column_to_rownames("barcode")

#mat <- mat[rowSums(mat) >= 5000, ]  # Keep only samples with at least 5000 reads


# assess variation in sampling effort, plotting sample effort curves
curve<- rarecurve(mat, step=1000, tidy=TRUE)

rare_soil<-Blast_ID%>%mutate(barcode=as.factor(barcode)) %>%
  right_join( curve%>%rename(barcode=Site,Samples=Sample),by = join_by('barcode'))%>% 
  ggplot(aes(x=Samples, y=Species, group=barcode,text=Site)) + 
  geom_line()+
  labs(tag = "(a)") +
 theme_classic()+
 theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18),
        plot.tag = element_text(size=25,hjust = 0.5),
        legend.position="none")+
  labs(y='Fungal OTUs', x= NULL) 

  
rare_soil
plotly::ggplotly(p, tooltip = "text")


Samples_included<-curve %>% 
  select(Site) %>%  distinct() %>% pull()


Guild_dat<-Guild_dat%>%  
  filter(barcode %in% Samples_included)

Decomp_Site<-Guild_dat%>%
  right_join(Blast_ID) %>% filter(!is.na(OTU))

Decomp_Site %>% distinct(barcode)

write.csv(Decomp_Site,'Soil_ITS_scripts/Jeff_proc/Processed_data/Updated_Soil_data.csv', row.names=FALSE)

write.csv(Blast_ID,'Soil_ITS_scripts/Jeff_proc/Processed_data/Site_Info.csv', row.names=FALSE)


Guild_dat_hp<-Guild_dat%>%
  filter(confidenceRanking=='Highly Probable')%>%
  left_join(Blast_ID)

write.csv(Guild_dat_hp,'Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_data_all_guild.csv', row.names=FALSE)

##########################
#Summary stats

# Load sample quality summary
Summary <- read_tsv('Processed_data/CG/ITS_data_quality_summary.tsv')

# Load summarised OTU table (likely aggregated by sample or taxonomy)
output <- read_tsv('Processed_data/CG/ITS_output_summarised.tsv')

# Filter out poorly sequenced samples (Sites 49 Tran 2 to Site 63)
description <- Summary %>%
  filter(barcode %in% Samples_included)

# Summarize read counts after rare OTU removal across remaining samples
description %>%
  summarise(
    total_reads = sum(n_after_removing_rare_otus),
    min_reads   = min(n_after_removing_rare_otus),
    max_reads   = max(n_after_removing_rare_otus)
  )

######################################

# Total read count in cleaned data (all samples), across all OTUs
total_reads<-clean_dat %>%
  filter(barcode %in% Samples_included) %>%
  summarise(total_reads = sum(count, na.rm = TRUE))

# Get min and max total read counts per sample
clean_dat %>%
  filter(barcode %in% Samples_included) %>%
  group_by(barcode) %>%
  summarise(total_count = sum(count)) %>%
  summarise(
    min_count = min(total_count),
    max_count = max(total_count)
  )

# Count distinct OTUs that have a genus assignment
t <- clean_dat %>%
  filter(barcode %in% Samples_included) %>%
  filter(!is.na(genus)) %>%
  distinct(OTU)

nrow(t)  # Returns number of distinct OTUs with genus assignment
 
# OPTIONAL: Count distinct OTUs with genus assignment *in the mycorrhizal dataset*
myco_dat %>%
  filter(barcode %in% Samples_included) %>%
  filter(!is.na(genus)) %>%
  distinct(OTU)

######################################
# 🔽 NEW STEP: Calculate total number of mycorrhizal reads in included samples

mycorrhizal_reads <- myco_dat %>%
  filter(barcode %in% Samples_included) %>%
  summarise(total_mycorrhizal_reads = sum(count, na.rm = TRUE))

mycorrhizal_reads

(mycorrhizal_reads/total_reads)*100
