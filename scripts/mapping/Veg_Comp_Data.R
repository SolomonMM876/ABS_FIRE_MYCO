library(readxl)
library(tidyr)
library(dplyr)
library(tidyverse)


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
  select(original_name,suggested_name,aligned_name,genus:taxon_rank)%>%
  filter(!is.na(suggested_name))

# mismatched_spp <- upd_spp_list %>%
#   filter(original_name != aligned_name, !is.na(aligned_name))

#######Austraits###########
#remotes::install_github("traitecoevo/austraits", dependencies = TRUE, upgrade = "ask")

#library(austraits)

#get_versions()[1,]

#austraits<-load_austraits(version = "6.0.0", path = "data/austraits")

# lookup_trait(austraits, "fire")
# lookup_trait(austraits, "root")



#resprouting_capacity#
#"N fixing"
# fire_recruit<-austraits%>%extract_trait("resprouting_capacity")
# 
# fire_recruit_trait<-flatten_database(fire_recruit)%>%select(taxon_name,genus,family,trait_name,value,taxon_rank,description, source_primary_citation)
# 
# fire_recruit_trait_clean<-fire_recruit_trait%>%
#   filter(trait_name=='resprouting_capacity')%>%
#   group_by(taxon_name) %>%
#   mutate(value = if_else(n_distinct(value) > 1, "resp_varies", first(value))) %>%
#   summarise(across(everything(), first), .groups = "drop")
# 
# upd_spp_list %>%rename(taxon_name=aligned_name)%>% select(-original_name,-suggested_name)%>%distinct()%>%group_by(taxon_name)%>%
#   mutate(n=n())->t
# 
# 
# species_fire_resp <- upd_spp_list %>%rename(taxon_name=aligned_name)%>% select(-original_name,-suggested_name)%>%distinct()%>%
#   left_join(
#     fire_recruit_trait_clean %>%
#       rename(fire_response = value))%>%
#   mutate(fire_response = case_when(
#     taxon_rank %in% c('genus', 'family') ~ 'insufficient_ID',  # Fixed typo in 'insufficient'
#     !is.na(fire_response) ~ fire_response,  # Keep existing values if they exist
#     TRUE ~ 'unknown'  # Assign 'unknown' if no other condition is met
#   ))
# 
# 
# 
# write.csv(species_fire_resp, 'Processed_data/Species_Fire_Resp.csv', row.names = FALSE)



#Clark paper
# 
# 
# clark_resprout<-read_excel('processed_data/clark_etal_2015.xlsx',skip=1)%>%rename(resprouter= `R +`, seeder= `S +`,taxon_name=`species records`)
# 
# fire_response<-clark_resprout%>%filter(Ecosystem=='Eucalypt Forest')%>%
#   mutate(taxon_name = str_replace(taxon_name, "_",' ')) %>%
#   group_by(taxon_name) %>%
#   distinct() %>%
#   mutate(n = n()) %>%
#   filter(!(n == 2 & resp.type == "X"))%>%
#   select(-n)
# 
# species_resprout <- upd_spp_list %>%rename(taxon_name=aligned_name)%>%
#   left_join(fire_response) 
# 
# 
# ABS_Floristics_fire_response<- ABS_Floristics %>%
#   filter(!ScientificName%in% c('Unknown','unknown','Unknown L'))%>%
#   left_join(fire_response, by =c("ScientificName"="taxon_name")) 
# 
# 
# fire_response_abun <- ABS_Floristics_fire_response %>%
#   mutate(
#     Transect1 = rowSums(select(., Q1.:Q15.), na.rm = TRUE),
#     Transect2 = rowSums(select(., Q16.:Q30.), na.rm = TRUE)
#   ) %>%
#   group_by(Site,resprouter) %>%
#   summarise(
#     Total_Transect1 = sum(Transect1),
#     Total_Transect2 = sum(Transect2))%>%
#   pivot_wider(names_from = resprouter, values_from = c(Total_Transect1, Total_Transect2), values_fill = 0) %>%
#   mutate(
#     resprouter_non_Transect1 = Total_Transect1_resprouter / (Total_Transect1_non_resprouter+Total_Transect1_resprouter+Total_Transect1_unsure),
#     resprouter_non_Transect2 = Total_Transect2_resprouter / (Total_Transect2_non_resprouter+Total_Transect2_resprouter+Total_Transect2_unsure))%>%
#   select(-starts_with('Total')) %>%
#   pivot_longer(cols = starts_with("resprouter"), 
#                names_to = c("resprouter_ratio", "Transect"),
#                names_pattern = "resprouter_to_non_(.*)(\\d)", 
#                values_to = "perc_resprouter_freq") %>%
#   select(-resprouter_ratio)
# 
# 
# 
# write.csv(species_resprout, file='Processed_data/resprout_seeder_abundance.csv',row.names=FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #"N fixing"
# N_fix<-austraits%>%extract_trait("nitrogen_fixing")
# 
# 
# N_fix_trait<-flatten_database(N_fix)%>%select(taxon_name,genus,family,trait_name,value,taxon_rank,description)
# 
# 
# 
# species_N_fix <- upd_spp_list %>%rename(taxon_name=aligned_name)%>%
#   left_join(
#     N_fix_trait %>%
#       rename(nitrogen_fixing = value) %>%
#       group_by(taxon_name) %>%
#       summarise(N_fixer_spp = paste(unique(nitrogen_fixing), collapse = " | "), .groups = "drop")    )
# 
# 
# 
# genus_N_fix<-species_N_fix%>%
#   left_join(N_fix_trait %>%
#               rename(nitrogen_fixing = value) %>%
#               group_by(genus) %>%
#               summarise(N_fixer_gen = paste(unique(nitrogen_fixing), collapse = " | "), .groups = "drop"))%>%
#   mutate(N_fixer= coalesce(N_fixer_spp,N_fixer_gen))%>% 
#   mutate(N_fixer= case_when(  grepl("\\|", N_fixer)|is.na(N_fixer)~ 'unsure',
#                               TRUE~N_fixer))
# 
# ABS_Floristics_N_fix <- ABS_Floristics %>%
#   filter(!ScientificName%in% c('Unknown','unknown','Unknown L'))%>%
#   left_join(genus_N_fix, by =c("ScientificName"="original_name")) 
#  
# 
# N_fixer_abun <- ABS_Floristics_N_fix %>%
#   mutate(
#     Transect1 = rowSums(select(., Q1.:Q15.), na.rm = TRUE),
#     Transect2 = rowSums(select(., Q16.:Q30.), na.rm = TRUE)
#   ) %>%
#   group_by(Site,N_fixer) %>%
#   summarise(
#     Total_Transect1 = sum(Transect1),
#     Total_Transect2 = sum(Transect2))%>%
#   pivot_wider(names_from = N_fixer, values_from = c(Total_Transect1, Total_Transect2), values_fill = 0) %>%
#   mutate(
#     N_fix_to_non_Transect1 = Total_Transect1_nitrogen_fixer / (Total_Transect1_non_nitrogen_fixer+Total_Transect1_nitrogen_fixer+Total_Transect1_unsure),
#     N_fix_to_non_Transect2 = Total_Transect2_nitrogen_fixer / (Total_Transect2_non_nitrogen_fixer+Total_Transect2_nitrogen_fixer+Total_Transect2_unsure))%>%
#   select(-starts_with('Total')) %>%
#   pivot_longer(cols = starts_with("N_fix"), 
#                names_to = c("N_fix_ratio", "Transect"),
#                names_pattern = "N_fix_to_non_(.*)(\\d)", 
#                values_to = "perc_N_fix_freq") %>%
#   select(-N_fix_ratio)
# 
# write.csv(N_fixer_abun, file='Processed_data/N_fixer_abundance.csv',row.names=FALSE)
# 
# 
# 
# 




######################################mycorr type#####
library(stringr)
brundrett_2017 <- read_csv("Processed_data/brundrett_2017_myco.csv")%>%
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
  select(genus,family,suggested_name,Myco_type)%>%distinct()



############################################
#calc frquency of mycorrhizal and non-mycorrhizal * and unknown

ABS_Floristics_root_str <- ABS_Floristics %>%
  filter(!ScientificName%in% c('Unknown','unknown'))%>%
  left_join(root_trait_tax, by =c("ScientificName"="original_name"))%>%
  filter(!Myco_type =='Unkown')%>%
  mutate(Myco_type = case_when(
    grepl("VAM|ECM|Ericoid|Orchid|MH|MHOrchid|MHVAM", Myco_type) ~ "Myco", 
    grepl("NM", Myco_type) ~ "Non_Myco", 
    TRUE ~ "Unknown",  # Covers NA or unspecified cases
    ))

# 

Myco_abun <- ABS_Floristics_root_str %>%
  mutate(
    Transect1 = rowSums(select(., Q1.:Q15.), na.rm = TRUE),
    Transect2 = rowSums(select(., Q16.:Q30.), na.rm = TRUE)
  ) %>%
  group_by(Site,Myco_type) %>%
  summarise(
    Total_Transect1 = sum(Transect1),
    Total_Transect2 = sum(Transect2))%>%
  pivot_wider(names_from = Myco_type, values_from = c(Total_Transect1, Total_Transect2), values_fill = 0) %>%
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
  group_by(Site, Myco_type) %>%
  summarise(
    Total_Transect1 = sum(Transect1),
    Total_Transect2 = sum(Transect2)
  ) %>%
  pivot_wider(names_from = Myco_type, values_from = c(Total_Transect1, Total_Transect2), values_fill = 0) %>%
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
###################################
#AM and ECM#########

# Reclassify Myco_type to separate VAM and ECM
ABS_Floristics_root_str_AM_ECM <- ABS_Floristics %>%
  filter(!ScientificName %in% c('Unknown', 'unknown')) %>%
  left_join(root_trait_tax, by = c("ScientificName" = "original_name")) %>%
  filter(!Myco_type == 'Unknown') %>%
  mutate(Myco_type = case_when(
    grepl("VAM", Myco_type) ~ "VAM",
    grepl("ECM", Myco_type) ~ "ECM",
    grepl("VAM|ECM|Ericoid|Orchid|MH|MHOrchid|MHVAM", Myco_type) ~ "Other_Myco",
    TRUE ~ "Non_Myco"
  ))

# Calculate total abundance per mycorrhizal type per site/transect
Myco_abun_AM_ECM <- ABS_Floristics_root_str_AM_ECM %>%
  mutate(
    Transect1 = rowSums(select(., Q1.:Q15.), na.rm = TRUE),
    Transect2 = rowSums(select(., Q16.:Q30.), na.rm = TRUE)
  ) %>%
  group_by(Site, Myco_type) %>%
  summarise(
    Total_Transect1 = sum(Transect1),
    Total_Transect2 = sum(Transect2),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Myco_type, values_from = c(Total_Transect1, Total_Transect2), values_fill = 0) %>%
  mutate(
    Prop_VAM_Transect1 = Total_Transect1_VAM / (Total_Transect1_VAM + Total_Transect1_ECM + Total_Transect1_Other_Myco+Total_Transect1_Non_Myco),
    Prop_VAM_Transect2 = Total_Transect2_VAM / (Total_Transect2_VAM + Total_Transect2_ECM + Total_Transect2_Other_Myco+ Total_Transect2_Non_Myco ),
    Prop_ECM_Transect1 = Total_Transect1_ECM / (Total_Transect1_VAM + Total_Transect1_ECM + Total_Transect1_Other_Myco + Total_Transect1_Non_Myco),
    Prop_ECM_Transect2 = Total_Transect2_ECM / (Total_Transect2_VAM + Total_Transect2_ECM + Total_Transect2_Other_Myco+ Total_Transect2_Non_Myco)
  ) %>%
  select(Site, starts_with("Prop_")) %>%
  pivot_longer(
    cols = starts_with("Prop_"),
    names_to = c("Myco_type", "Transect"),
    names_pattern = "Prop_(VAM|ECM)_Transect(.*)",
    values_to = "perc_myco_type_freq"
  )%>%
  mutate(Transect=str_remove(Transect,'Transect'))%>%
  pivot_wider(names_from = Myco_type, values_from = perc_myco_type_freq, values_fill = 0)

# Join with total mycorrhizal host abundance
Myco_plant_spp <- Myco_plant_spp%>%
  left_join(Myco_abun_AM_ECM, by = c("Site", "Transect"))

write.csv(Myco_plant_spp, file='Processed_data/Myco_host_abundance.csv',row.names=FALSE)

library(ggplot2)
library(ggpubr)

# Scatter plot with regression line for VAM_to_myco
plot_vam <- ggplot(Myco_plant_spp, aes(x = VAM_to_myco, y = perc_myco_host_freq)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "blue") +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +
  labs(title = "Relationship between VAM to Myco and Mycorrhizal Frequency",
       x = "VAM to Myco Ratio", y = "Percentage Mycorrhizal Hosts") +
  theme_minimal()

# Scatter plot with regression line for ECM_to_myco
plot_ecm <- ggplot(Myco_plant_spp, aes(x = ECM_to_myco, y = perc_myco_host_freq)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red") +
  stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~"))) +
  labs(title = "Relationship between ECM to Myco and Mycorrhizal Frequency",
       x = "ECM to Myco Ratio", y = "Percentage Mycorrhizal Hosts") +
  theme_minimal()

# Print plots
print(plot_vam)
print(plot_ecm)

# Perform linear regression and check significance
lm_vam <- lm(perc_myco_host_freq ~ VAM_to_myco, data = Myco_plant_spp)
lm_ecm <- lm(perc_myco_host_freq ~ ECM_to_myco, data = Myco_plant_spp)

# Summarize regression results
summary(lm_vam)  # R² and p-value for VAM_to_myco
summary(lm_ecm)  # R² and p-value for ECM_to_myco
