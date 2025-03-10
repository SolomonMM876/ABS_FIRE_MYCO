library(readxl)
library(tidyr)
library(dplyr)


# Load and process vegetation and site information data
VEG_Transect <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                           sheet = "Transect.Level_Data") %>%
  mutate(
    Site = gsub('ABS00|ABS0', "", Site),
    Transect = gsub('T', "", Transect)
  )

Study_sites<-read.csv('Processed_data/All_Bag_data.csv')%>%
  select(Site)%>%distinct()%>%mutate(Site=as.character(Site))


ABS_Floristics<- read.csv("Raw_data/Site_Data/ABS_Floristic plot data_ANALYSES_V2_20240819_CEG.cleaned.csv")[,-1]%>%
  mutate(SiteNo = gsub('ABS00|ABS0', "", SiteNo))%>%
  rename(Site=SiteNo)#%>%
  #inner_join(Study_sites, by='Site')

library(APCalign)
#load for austraits
tax_resources <- load_taxonomic_resources()

spp_list<-ABS_Floristics%>%distinct(ScientificName)

spp_list<-align_taxa(spp_list$ScientificName, resources = tax_resources)

upd_spp_list<-update_taxonomy(spp_list, taxonomic_splits = "most_likely_species",resources = tax_resources)%>%
  rename(taxon_name=`suggested_name`)%>%
  select(original_name,taxon_name,genus:taxon_rank)%>%
  filter(!is.na(taxon_name))

# mismatched_spp <- upd_spp_list %>%
#   filter(original_name != aligned_name, !is.na(aligned_name))

#######Austraits, not using this###########
# library(austraits) 
# 
# get_versions()[1,]
# 
# austraits<-load_austraits(version = "6.0.0", path = "data/austraits")
# 
# #lookup_trait(austraits, "fire")
# #lookup_trait(austraits, "root")
# 
# #"post_fire_recruitment" 
# root_str<-austraits%>%extract_trait("root_structure")
# 
# 
# #first we summarise by taxon
# Root_trait<-flatten_database(root_str)%>%select(taxon_name,genus,family,trait_name,value,taxon_rank)
# 
######################################
library(stringr)
brundrett_2017 <- read_csv("Processed_data/brundrett_20170327.csv")%>%
  #correct family names so joins work
  mutate(Family = str_replace(Family, "^(Euphorbaceae|Fabaceae)\\b.*", "\\1"),
        Family = str_replace(Family, "Euphorbaceae", "Euphorbiaceae"))

 
# First join: Match trait data based on Genus, some taxa are dual, but we just care if they are or are not
genus_myco <- upd_spp_list %>%
  left_join(
    brundrett_2017 %>%
      rename(Myco_type = Primary.Obs) %>%
      group_by(Genus) %>%
      summarise(Myco_type = paste(unique(Myco_type), collapse = " | "), .groups = "drop"),  
    by = c("genus" = "Genus")
  ) %>%
  mutate(
    Myco_type = case_when(
      genus %in% c("Thysanotus", "Cassytha", "Olax") ~ "NM",
      genus %in% c('Doryanthes',"Brunoniella",'Caesia') ~ "Unkown",
      TRUE ~ Myco_type),
    Notes = case_when(
      genus == "Thysanotus" ~ "Unique sub-epidermal association - McGee 1988, Brundrett & Abbott 1991.",
      genus == "Cassytha" ~ "Can have members of family that are mycorrhizal, but all species observed are hemi-parasitic growing on stems and are NM as such.",
      genus == "Brunoniella" ~ "Type not determined due to lack of records, removed from analysis.",
      #taxon_name == "Doryanthes excelsa" ~ "Type not determined due to lack of records, classified as NM for analysis.",
      #taxon_name == "Caesia sp." ~ "Type not determined due to lack of records, classified as NM for analysis.",
      genus == "Olax" ~ "Non-mycorrhizal Bellgard et al. 1994",
      is.na(Myco_type) ~ "Classification based on other genera in family",
      TRUE ~ ""))

family_myco<-genus_myco%>%
  left_join(brundrett_2017%>% rename(Myco_type_fam=Primary.Obs)%>%
               group_by(Family) %>%  
               summarise(Myco_type_fam = paste(unique(Myco_type_fam), collapse = " | "), .groups = "drop",),
            by= c('family'='Family'))%>%
  mutate(Myco_type= coalesce(Myco_type,Myco_type_fam))%>%
  select(-Myco_type_fam)

root_trait_tax<-family_myco%>%
  left_join(genus_myco)

myco_type_pub<-root_trait_tax%>%rename(Species=original_name, Family=family,Genus=genus,Type=Myco_type)%>%
  select(Genus,Family,Type,Notes)%>%distinct()

write.csv(myco_type_pub,'processed_data/Myco_typ_pub.csv', row.names = FALSE)

unk_myco_plants<-root_trait_tax%>%
  filter(is.na(Myco_type))%>%
  select(genus,family,taxon_name,Myco_type)%>%distinct()



############################################
#calc frquency of mycorrhizal and non-mycorrhizal * and unknown

ABS_Floristics_root_str <- ABS_Floristics %>%
  filter(!ScientificName%in% c('Unknown','unknown'))%>%
  left_join(root_trait_tax, by =c("ScientificName"="original_name"))%>%
  filter(Myco_Status=='Unkown')%>%
  mutate(Myco_Status = case_when(
    grepl("VAM|ECM|Ericoid|Orchid|MH|MHOrchid|MHVAM", Myco_type) ~ "Myco", 
    grepl("NM", Myco_type) ~ "Non_Myco", 
    TRUE ~ "Unknown"  # Covers NA or unspecified cases
  ))

# 

Myco_abun <- ABS_Floristics_root_str %>%
  mutate(
    Transect1 = rowSums(select(., Q1.:Q15.), na.rm = TRUE),
    Transect2 = rowSums(select(., Q16.:Q30.), na.rm = TRUE)
  ) %>%
  group_by(Site,Myco_Status) %>%
  summarise(
    Total_Transect1 = sum(Transect1),
    Total_Transect2 = sum(Transect2))%>%
  pivot_wider(names_from = Myco_Status, values_from = c(Total_Transect1, Total_Transect2), values_fill = 0) %>%
  mutate(
    Myco_to_non_Transect1 = Total_Transect1_Myco / (Total_Transect1_Non_Myco+Total_Transect1_Myco),
    Myco_to_non_Transect2 = Total_Transect2_Myco / (Total_Transect2_Non_Myco+Total_Transect2_Myco))%>%
  select(-starts_with('Total')) %>%
  pivot_longer(cols = starts_with("Myco_to"), 
               names_to = c("Myco_ratio", "Transect"),
               names_pattern = "Myco_to_non_(.*)(\\d)", 
               values_to = "perc_myco_host_freq") %>%
  select(-Myco_ratio)

Myco__total <- ABS_Floristics_root_str %>%
  mutate(
    Transect1 = rowSums(select(., Q1.:Q15.), na.rm = TRUE),
    Transect2 = rowSums(select(., Q16.:Q30.), na.rm = TRUE)
  ) %>%
  group_by(Site, Myco_Status) %>%
  summarise(
    Total_Transect1 = sum(Transect1),
    Total_Transect2 = sum(Transect2)
  ) %>%
  pivot_wider(names_from = Myco_Status, values_from = c(Total_Transect1, Total_Transect2), values_fill = 0) %>%
  select(-ends_with("Non_Myco"),-ends_with("Unknown")) %>%
  pivot_longer(
    cols = starts_with("Total_"),
    names_to = "Transect",
    names_pattern = "Total_Transect(\\d+)_Myco",
    values_to = "myco_host_total"
  ) 

Myco_plant_spp<-Myco__total%>%
  left_join(Myco_abun)





cor(Myco_plant_spp$myco_host_total,Myco_plant_spp$perc_myco_host_freq)









write.csv(Myco_plant_spp, file='Processed_data/Myco_host_abundance.csv',row.names=FALSE)

