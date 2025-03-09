library(tidyverse)
library(dplyr)
library(readr)

dat <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_clean.tsv')

#filter out data for CNP measurements and community
CNP_dat<-dat%>%
  filter(is.na(Site))

write.csv(CNP_dat, file='Processed_data/CNP_seq_dat.csv',row.names=FALSE)

#extract all unknown taxa from sequencing
tax<-dat%>%
  filter(is.na(guild)) %>%
  select(OTU,kingdom:guild)%>%
  distinct()


write.csv(tax, file='Processed_data/unknown_taxa_bags.csv',row.names=FALSE)


clean_dat <- dat %>%
  #remove transect pooled samples for CN
  anti_join(CNP_dat, by = colnames(CNP_dat)) %>%
  #Samples 95 and 96 were pooled. I have removed all separated them,
  #but for transect analyses, but they should be removed when looking at location based correlations
  mutate(Tube_ID = case_when(
    Location == 3 & Transect == 2 & Site == 11 ~ "95",
    Location == 33 & Transect == 2 & Site == 11 ~ "96",
    TRUE ~ Tube_ID  # Keep the original value if no condition is met
  ))%>%
  filter(!Tube_ID %in% c('84'))#Tube 84 is linked to two locations 11-1-16 and 11-2-16 due to extraction error
 
dat%>%
  group_by(barcode,Tube_ID)%>%
  summarise()->temp

Guild_dat<-clean_dat%>%
  mutate(guild2 = case_when(
    str_detect(guild, "Saprotroph") ~ "sap",
    guild == "Arbuscular Mycorrhizal" ~ "amf",
    guild == "Ectomycorrhizal" ~ "ecm",
    guild ==  "Plant Pathogen" ~ "path",
    guild == 'Endophyte'~ "endo",
    TRUE ~ 'other_unknown'   )) %>%  # Assigns NA to any rows that do not match the specified conditions
  group_by(Tube_ID, guild2) %>% 
  summarise(count=sum(count))%>%
  mutate(Tube_ID=as.numeric(Tube_ID))


library(fungaltraits)

tax<-clean_dat%>%
  select(OTU,kingdom:guild)%>%
  distinct()

fungal_traits<-fungal_traits()

traits_to_join <- c("melanin_content", "guild_fg", "em_text", 
                    "em_expl","spore_size","melanin_content")
# Join by Species
tax_species_traits <- tax %>%
  left_join(fungal_traits %>%
              select(species, all_of(traits_to_join))%>%distinct(), by = "species")

# Join by Genus
clean_dat_genus <- tax %>%
  left_join(fungal_traits
          %>% select(Genus, all_of(traits_to_join))%>%distinct(), 
    by = c("genus" = "Genus"), 
    suffix = c("_species", "_genus"))





Bag_data<-read.csv('Processed_data/All_Bag_Site_Info.csv')%>%
  dplyr::select(Tube_ID,Site,Transect,Fire.Severity,Fire.Interval, Location,myc_2nd_w_est_yield,Days_Installed)


Bag_Guild<-Bag_data%>%right_join(Guild_dat, by="Tube_ID" )




Bag_Guild %>%
  group_by(site_transect = paste("Site", Site, "T", Transect)) %>%
  mutate(
    Site = as.factor(Site),
    guild2 = factor(guild2, levels = c(setdiff(unique(guild2), "other_unknown"), "other_unknown")),
    # Calculate relative abundance within each site/transect
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
    fill = "Guild"
  ) +
  scale_y_continuous(limits = c(0, 100))
  
  
  #code for biomass calcs TO BE DELT WITH LATER
  # mutate(#mass extraacted????????
          #sap_biomass_bag= (myc_2nd_w_est_yield*ra_sap),
  #       seq_corrected_myc= myc_2nd_w_est_yield-sap_biomass_bag,
  #       seq_corrected_myc= if_else(seq_corrected_myc==0, 0.00114/2, seq_corrected_myc),#lowest non-zero value in biomass produced
  #        seq_cor_Biomass_day= seq_corrected_myc/ Days_Installed,
  #        log10_seq_cor_myc =log10(seq_corrected_myc),
  #        log10_cor_biomass_day = log10(seq_cor_Biomass_day),
  #        cor_biomass_g_ha_day = seq_cor_Biomass_day * (1e+06 / 15),  # Convert to g/ha/day
  #          )


#write_xlsx(Bag_Seq, path='Processed_data/corrected_biomass_df.xlsx')
write.csv(clean_dat, file='Processed_data/cleaned_seq_dat.csv',row.names=FALSE)


dat_summary <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_summarised.tsv')


