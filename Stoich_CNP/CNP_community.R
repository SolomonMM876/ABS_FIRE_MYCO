library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(stringr)
library(readr)


Bag_data<-read.csv('Processed_data/All_Bag_Site_Info.csv')%>%
  dplyr::select(Tube_ID,Site,Transect,Fire.Severity,Fire.Interval)%>%
  group_by(Site,Transect)

CNP_Seq<-read.csv('Processed_data/CNP_seq_dat.csv') %>%
  mutate(Tube_ID = if_else(Tube_ID == "26.1A26.2", "26.1", Tube_ID))%>%
  #26.2 was added to 26.1, but contamination was minimal therefore 26.1A26.2 is transformed to just 26.1
  separate(Tube_ID, into = c("Site", "Transect"), sep = "\\.")%>%
  mutate(Site=as.numeric(Site),
         Transect=as.numeric(Transect))

dat_summary_Jeff <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_summarised.tsv')%>%
  filter(str_detect(barcode, "T"))

Jeff_alpha<-Bag_data%>%select(-Tube_ID)%>%distinct()%>%
  right_join(dat_summary_Jeff%>%left_join(CNP_Seq%>%select(barcode,Site,Transect)%>%distinct()))%>%select(-sample)

#calc guild proportion per sample
CNP_Guild<-CNP_Seq%>%
  mutate(guild2 = case_when(
    str_detect(guild, "Saprotroph") ~ "sap",
    guild == "Arbuscular Mycorrhizal" ~ "amf",
    guild == "Ectomycorrhizal" ~ "ecm",
    guild == "Plant Pathogen" ~ "path",
    guild == 'Endophyte' ~ "endo",
    guild == 'Orchid Mycorrhizal' ~ 'orch_mf',
    guild == 'Ericoid Mycorrhizal' ~ 'Ercm',
    is.na(guild) ~ 'no_guild',  # Assigns 'no_guild' if guild is NA
    TRUE ~ 'other_unknown'      # Handles other non-matching values
  )) %>%  # Assigns NA to any rows that do not match the specified conditions
  group_by(Site,Transect,guild2) %>% 
  summarise(count=sum(count))%>%
  ungroup() %>%  # Remove grouping for easier calculations
  # Calculate total reads per sample (Site, Transect)
  group_by(Site, Transect) %>%
  mutate(total_reads = sum(count)) %>%
  ungroup() %>%
  # Calculate percentage of total reads per guild
  mutate(percentage = (count / total_reads) * 100,
         Group = case_when(
           guild2 %in% c("amf", "ecm", "orch_mf", "Ercm") ~ "Mycorrhizal",
           guild2 %in% c("sap","path", "other_unknown")  ~ 'Non-Mycorrhizal',
           TRUE ~ "Unknown"
         ))

Bag_comp<-CNP_Guild%>%group_by(Group)%>%summarise()


Sample_Guild<-Bag_data%>%right_join(CNP_Guild, by = c('Site','Transect'),relationship = "many-to-many")


Sample_Guild %>%
  group_by(Sample = paste("Site", Site, "T", Transect)) %>%
  mutate(
    Site = as.factor(Site),
    guild2 = factor(guild2, levels = c(setdiff(unique(guild2), "other_unknown"), "other_unknown")),
    # Calculate relative abundance within each site/transect
    rel_abundance = count/sum(count) * 100,
    # Create a factor with custom ordering for site/transect combinations
    Sample = factor(Sample,
                           levels = unique(Sample)[order(Site, Transect)]))%>%
  ggplot(aes(x = Sample, y = rel_abundance, fill = guild2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Fire.Interval, scales = "free_x") +
  labs(
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = "Guild"
  ) +
  scale_y_continuous(limits = c(0, 100))






##### Prep data



#calc guild proportion per sample
CNP_myco_comm<-CNP_Seq%>%
  mutate(guild2 = case_when(
    str_detect(guild, "Saprotroph") ~ "sap",
    guild == "Arbuscular Mycorrhizal" ~ "amf",
    guild == "Ectomycorrhizal" ~ "ecm",
    guild == "Plant Pathogen" ~ "path",
    guild == 'Endophyte' ~ "endo",
    guild == 'Orchid Mycorrhizal' ~ 'orch_mf',
    guild == 'Ericoid Mycorrhizal' ~ 'Ercm',
    is.na(guild) ~ 'no_guild',  # Assigns 'no_guild' if guild is NA
    TRUE ~ 'other_unknown'      # Handles other non-matching values
  ))  %>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))%>%
  filter(guild %in%c("Ectomycorrhizal","Arbuscular Mycorrhizal","Orchid Mycorrhizal","Ericoid Mycorrhizal"))%>%
  mutate(barcode=str_remove(barcode,'A2'))# fix note about DNA extraction
  
write.csv(CNP_myco_comm, 'Processed_data/CNP_seq__myco_dat.csv', row.names = FALSE)