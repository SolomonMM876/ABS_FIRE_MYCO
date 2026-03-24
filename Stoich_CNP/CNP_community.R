library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(stringr)
library(readr)

#source('df_joins.R')

# 
# Bag_data <- read.csv('Processed_data/All_Bag_data.csv') %>%
#   mutate(Tube_ID = as.factor(paste(Site, Transect, sep = ".")),
#          Tube_ID=as.factor(Tube_ID))%>%
#   group_by(Tube_ID,Site, Transect,Fire.Severity,Fire.Interval) %>%
#   summarize(across(
#     c(log10_biomass_day, pH, Ortho_P_mg_kg, Ammonia_mg_kg, Nitrate_mg_kg, 
#       #Sample_mg_CN, C_N, Sample_mg_P, C_P, N_P, Carb_Hyph,Nitrog_Hyph,Phos_Hyph,
#       Longitude,Latitude,
#       myco_host_total, perc_myco_host_freq, perc_N_fix_freq), 
#     ~ mean(.x, na.rm = TRUE)
#   )) %>%
#   distinct()

CNP_clean<-read.csv("Processed_data/CNP_clean.csv") %>% 
  rename(Tube_ID=Sample) %>% 
  mutate(Tube_ID=as.factor(Tube_ID))



CNP_Seq<-read.csv('Processed_data/CNP_seq_dat.csv') 

total_reads_per_sample <- CNP_Seq %>%
  group_by(Tube_ID) %>%
  summarize(total_reads = sum(count, na.rm = TRUE))%>%
  mutate(Tube_ID=as.factor(Tube_ID))

#calc guild proportion per sample
CNP_Guild<-CNP_Seq%>%select(-c(sample,Site,Transect,Location,Note))%>%mutate(Tube_ID=as.factor(Tube_ID))%>%
  left_join(CNP_clean)%>%
  filter(confidenceRanking == "Highly Probable" | OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073")) %>%  # Keep "Highly Probable" & specified OTUs
  filter(!str_detect(guild, "-") | OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073")) %>%  # Exclude multiple guilds except for specified OTUs
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
                            TRUE ~ 'unassigned'),
         genus = case_when(  # Separate case_when() for genus assignment
           OTU == 'ITSall_OTUa_10308' ~ 'Tomentella',
           OTU == 'ITSall_OTUk_1642' ~ 'Coltricia',
           OTU == 'ITSall_OTUm_3073' ~ 'Ionosporus',
           TRUE ~ genus  # Keep the existing genus if OTU is not listed
         ),
         note = if_else(OTU %in% c("ITSall_OTUa_10308", "ITSall_OTUk_1642", "ITSall_OTUm_3073"), 
                        "Included due to manual BLAST confirmation", 
                        NA_character_)) %>% # Add note for manually included OTUs
  left_join(total_reads_per_sample, by = "Tube_ID") %>%  # Join total reads per sample
  group_by(Tube_ID) %>% 
  mutate(rel_abundance = (count / total_reads) * 100) %>%   # Calculate relative abundance based on all reads
  ungroup()

all_tax<-CNP_Guild %>% 
  select(OTU,kingdom:species,confidenceRanking, guild,guild2)%>%
  distinct()

temp<-CNP_Guild%>%
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
           guild2 %in% c("Arbuscular Mycorrhizal", "Ectomycorrhizal", "Orchid Mycorrhizal") ~ "Mycorrhizal",
           guild2 %in% c("Saprotroph")  ~ 'Saprotroph',
           guild2 %in% c("Endophyte")  ~ 'Endophyte'
         ))



CNP_Guild_RA<-CNP_Guild%>%
  group_by(Site,Transect,barcode)%>%
  summarise(reads_samp=sum(count))%>%
  left_join(CNP_Guild)%>%
  ungroup()%>%
  left_join(CNP_Guild%>%
              group_by(Fire.Severity)%>%
              summarise(severity_reads=sum(count)))%>%
  ungroup()%>%
  left_join(CNP_Guild%>%
              group_by(Fire.Interval)%>%
              summarise(interval_reads=sum(count)))%>%
  ungroup()%>%
  mutate(total_reads=sum(count),
         RA_samp= count/reads_samp,
         RA_total_reads= count/total_reads,
         RA_total_interval= count/interval_reads,
         RA_total_severity= count/severity_reads)




CNP_Guild%>%
  ggplot(aes(x = Tube_ID, y = rel_abundance, fill = guild2)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~Fire.Interval, scales = "free_x") +
  labs(
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = "Guild"
  ) 

# write.csv(CNP_Guild_RA, 'Processed_data/CNP_seq__myco_dat.csv', row.names = FALSE)
# 
# otu_table<-CNP_Guild%>%select(Tube_ID, OTU, count) %>% 
#   pivot_wider(names_from = OTU, values_from = count, values_fill = 0)%>%
#   column_to_rownames('Tube_ID')
# 
# 
# write.csv(otu_table,'Processed_data/CNP_otu_table_bag.csv')



#format required for some vegan functions
wide_myco<-CNP_Guild  %>%
  pivot_wider(
    names_from = OTU,         # Column containing OTU names that will become new columns
    values_from = count,      # Column containing values that will fill the new columns
    id_cols = Tube_ID,  # Column(s) to keep as identifier
    values_fill = 0
  )%>% mutate(Tube_ID=as.factor(Tube_ID))

Bag_Seq_wide<-CNP_clean%>% 
  left_join(CNP_Guild%>% 
            group_by(Tube_ID) %>% 
           summarise(total_count=sum(count)))%>%
  left_join(wide_myco %>% mutate(Tube_ID=as.factor(Tube_ID)))%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))%>%ungroup()

#write.csv(Bag_Seq_wide, 'Processed_data/CNP_Bag_Seq_wide.csv', row.names = FALSE)#bag data with selected meta and Myco communities
#write.csv(Bag_data, 'Processed_data/CNP_Updated_Bag_data.csv', row.names = FALSE)#Just Site info, veg nute meta etc
#write.csv(all_tax, 'Processed_data/CNP_Bag_dat_tax.csv', row.names = FALSE)#mycorrhizal taxa included





# next analysis - permanova

# first remove samples that have no  reads because they cause errors below
Bag_Seq_wide<-Bag_Seq_wide%>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0)) %>% 
  filter(!is.na(C_P))

# extract the community table, save as a new object
mat_myco<- Bag_Seq_wide %>% select(starts_with("ITSall"))

colnames(Bag_Seq_wide)

rda<-rda(mat_myco~ Carbon + Nitrogen + Hydrogen+ C_N+  Percent_Phos_ + C_P+ N_P + Condition (Site)
        , data=Bag_Seq_wide)

rda<-rda(mat_myco~ + Condition (Site)
         , data=Bag_Seq_wide)


plot(rda)









#Testing for differences in between samples
adonis2(mat_myco~ Site, data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)

anosim(vegdist(mat_myco, method = "bray"), Bag_Seq_wide$Site)

mrpp(vegdist(mat_myco, method = "bray"), Bag_Seq_wide$Site)


dispersion <- betadisper(vegdist(mat_myco, method = "robust.aitchison"), Bag_Seq_wide$Fire.Interval)
anova(dispersion)

cap.all <- capscale(mat_myco~ Condition (Site)
                    , data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)

cap.all <- capscale(mat_myco~ Carbon + Nitrogen + Hydrogen+ C_N+  Percent_Phos_ + C_P+ N_P +Condition (Site)
                    , data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)
anova(cap.all)
summary(cap.all)
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
cap.all

plot(cap.all)
proportions<-round(cap.all$CCA$eig/cap.all$tot.chi *100, 1) # proportion of variation associated with each axis


# produce a nice plot
# Extract scores
scrs <- scores(cap.all, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site<- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')
scrs_biplot <- scrs %>% filter(score=='biplot')

interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_colors <- c("Hgih" = "blue", "Low" = "green")

library(ggrepel)

# Custom color mappings
interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_shapes <- c("High" = 19, "Low" = 1)

# Combine site scores with metadata
site_data <- cbind(Bag_Seq_wide, scrs_site)

# Make the plot
p2 <- ggplot(site_data, aes(x = CAP1, y = CAP2)) + 
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +  
  geom_point(aes(colour = Fire.Interval, shape = Fire.Severity), size = 6, stroke = 1.2) + 
  geom_text(aes(label = Site), color = 'black', size = 3, vjust = -1) +
  scale_shape_manual(values = severity_shapes) +
  scale_colour_manual(values = interval_colors) +
  labs(
    x = paste0("CAP1 (", proportions[1], "%)"), 
    y = paste0("CAP2 (", proportions[2], "%)"), 
    colour = "Fire Frequency", 
    shape = "Fire Severity"
  ) +
  geom_segment(data = scrs_biplot %>% filter(abs(CAP1) > 0.3 | abs(CAP2) > 0.3),
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(type = "closed", length = unit(3, 'mm')),
               color = "black") +
  geom_text_repel(data = scrs_biplot %>% filter(abs(CAP1) > 0.3 | abs(CAP2) > 0.3),
                  aes(x = CAP1, y = CAP2, label = label),
                  size = 5, color = "black") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

p2
