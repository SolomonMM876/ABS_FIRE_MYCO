
source('Soil_ITS_Scripts/Jeff_Proc/Soil_ITS_data_prep_Jeff_up.R')
library(ggplot2)
library(patchwork)



dat_myco_RA_soil %>%
  mutate(
    Fire_Int_Sev = interaction(Fire.Interval, Fire.Severity, sep = " x "))

#this script is only for non interactive effects#############
# top_by_severity_soil <- dat_myco_RA_soil %>%
#   group_by(Fire.Severity, OTU) %>%
#   summarise(total_reads = sum(count, na.rm = TRUE),
#             percent=(total_reads/first(severity_reads))*100,
#             .groups = "drop") %>%
#   group_by(Fire.Severity) %>%
#   slice_max(order_by = percent, n = 5) %>%
#   arrange(Fire.Severity, desc(percent))%>% 
#   left_join(myco_tax %>% select(OTU, genus,species)) %>% 
#   mutate(source='soil')
# 
# top_by_severity_soil
# 
# 
# top_by_interval_soil<- dat_myco_RA_soil %>%
#   group_by(Fire.Interval, OTU) %>%
#   summarise(total_reads = sum(count, na.rm = TRUE),
#             percent=(total_reads/first(interval_reads))*100,
#             .groups = "drop") %>%
#   group_by(Fire.Interval) %>%
#   slice_max(order_by = percent, n = 5) %>%
#   arrange(Fire.Interval, desc(percent))%>% 
#   left_join(myco_tax %>% select(OTU, genus,species)) %>% 
#   mutate(source='soil')
# 
# 
# top_by_interval_soil
###############


#this is for looking at the interactive effects
# Step 1: Create interaction variable
soil_myco <- myco_dat %>%
  mutate(Fire_Int_Sev = interaction(Fire.Interval, Fire.Severity, sep = " x "))

# for both dfs
dat_myco_RA_soil <- dat_myco_RA_soil %>%
  mutate(Fire_Int_Sev = interaction(Fire.Interval, Fire.Severity, sep = " x "))


#join them together
dat_myco_RA_soil<-dat_myco_RA_soil%>%
  left_join(soil_myco%>%
              group_by(Fire_Int_Sev)%>%
              summarise(interaction_reads=sum(count)))


# Step 2: Summarise and get top OTUs by interaction
top_by_interaction_soil <- dat_myco_RA_soil %>%
  group_by(Fire_Int_Sev, OTU) %>%
  summarise(
    total_reads = sum(count, na.rm = TRUE),
    percent = (total_reads / first(interaction_reads)) * 100,
    .groups = "drop"
  ) %>%
  group_by(Fire_Int_Sev) %>%
  slice_max(order_by = percent, n = 5) %>%
  arrange(Fire_Int_Sev, desc(percent)) %>%
  left_join(myco_tax %>% select(OTU, genus, species), by = "OTU") %>%
  mutate(source = "soil")

# View result
top_by_interaction_soil





# Identify the top 10 most abundant genera
top_genera_soil <- dat_myco_RA_soil %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(RA_total_severity, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(genus)

# Create a genus_grouped factor with proper ordering
genus_order <- dat_myco_RA_soil %>%
  mutate(genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus)) %>%
  group_by(genus_grouped) %>%
  summarise(total_abundance = sum(RA_total_severity, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  pull(genus_grouped)

genus_order <- c(setdiff(genus_order, "Other"), "Other")

# Add interaction term for plotting
dat_myco_RA_soil <- dat_myco_RA_soil %>%
  mutate(
    genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus),
    genus_grouped = factor(genus_grouped, levels = genus_order),
    Fire_Int_Sev = interaction(Fire.Interval, Fire.Severity, sep = " x ")
  )


covariate_labels <- c(
  "Long x High"   = "Long\nx\nHigh",
  "Short x Low"    = "Short\nx\nLow",
  "Short x High"  =  'Short\nx\nHigh',
  "Long x Low"   =   'Long\nx\nLow')
  
  dat_myco_RA_soil %>% distinct(Fire_Int_Sev)


# Plot
Interaction_Genus_myco_soil <- dat_myco_RA_soil %>%
  group_by(Fire_Int_Sev) %>%
  summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>% 
  left_join(dat_myco_RA_soil) %>% 
  ggplot(aes(x = Fire_Int_Sev, y = count/total, fill = genus_grouped)) +
  geom_bar(stat = 'identity', position = position_stack(), width = 0.6) +
  scale_fill_viridis_d(option = "H", name = "Genus",direction = -1) +
  theme_classic()+
  scale_x_discrete(labels = covariate_labels) +
  labs(
    y = '',
    x = '',
    tag = "(a)"
  ) +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 25,lineheight = 0.7),
    axis.text.y = element_text(hjust = 0.5, size = 25),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 28),
    plot.tag = element_text(size = 22, hjust = 0.5),
    legend.position = "right"
  )

Interaction_Genus_myco_soil



library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)

fire_levels <- c(
  "Long x Low",
  "Long x High",
  "Short x Low",
  "Short x High"
)

top_by_interaction_soil %>% 
  distinct(Fire_Int_Sev)

# Step 1: Create label for plotting
top_soil <- top_by_interaction_soil %>% 
  mutate(
    OTU_clean = str_remove(OTU, "ITSall_"),
    OTU_genus = paste0(genus, " spp. ", OTU_clean),
    OTU_genus = fct_reorder(OTU_genus, percent)) %>% 
  mutate( Fire_Int_Sev = factor(Fire_Int_Sev, levels = fire_levels)) 

library(tidytext)


TOP_OTU_soil<-top_soil %>% 
  mutate(OTU_genus = reorder_within(OTU_genus, percent, Fire_Int_Sev)) %>%
  ggplot(aes(x = OTU_genus, y = percent)) +
  geom_col(aes(fill = OTU), width = 0.7) +
  coord_flip() +
  scale_x_reordered() +
  facet_wrap(~Fire_Int_Sev, ncol = 1, scales = "free_y") +
  scale_fill_viridis_d(option = "H", name = "Genus",direction = -1) +
  guides(fill = "none") +  # Remove legend
  scale_y_continuous(expand = c(.01,0)) +  # <-- removes space between bars and axis
  labs(
    x = '',
    y = "",
    tag= '(a)'
  ) +
  theme_classic() +
  theme(
    legend.position = 'top',
    axis.text.x = element_text(hjust = 0.5, size = 25),
    axis.text.y = element_text(hjust = 0.5, size = 25),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    strip.text = element_text(size = 20),
    plot.tag = element_text(size = 22, hjust = 0.5),
    panel.spacing = unit(1, "lines")
  )

TOP_OTU_soil





#######################################


library(dplyr)
library(ggplot2)

# # Step 1: Create interaction variable
# dat_myco_RA_soil <- dat_myco_RA_soil %>%
#   mutate(Fire_Int_Sev = interaction(Fire.Interval, Fire.Severity, sep = " x "))

# Step 2: Count frequency of OTUs (number of samples where OTU is present)

# Count total sites per Fire_Int_Sev
sample_totals <- dat_myco_RA_soil %>%
  distinct(Site,Transect, Fire_Int_Sev) %>%
  count(Fire_Int_Sev, name = "total_samps")

# Count OTU presence per site and join with site totals
top_by_sample_interaction_percent <- dat_myco_RA_soil %>%
  filter(count > 0) %>%
  group_by(Fire_Int_Sev, OTU) %>%
  summarise(sample_count = n_distinct(Site,Transect), .groups = "drop") %>%
  left_join(sample_totals, by = "Fire_Int_Sev") %>%
  mutate(percent = (sample_count / total_samps) * 100) %>%
  group_by(Fire_Int_Sev) %>%
  slice_max(order_by = percent, n = 5) %>%
  arrange(Fire_Int_Sev, desc(percent)) %>%
  left_join(myco_tax %>% select(OTU, genus, family, order, species, class), by = "OTU") %>%
  mutate(source = "soil") %>%
  ungroup()



# Step 3: Plot
top_freq_soil<-top_by_sample_interaction_percent %>%
  mutate(
    OTU_clean = str_remove(OTU, "ITSall_"),
    ID = case_when(
      !is.na(genus) ~ genus,
      !is.na(family) ~ family,
      !is.na(order) ~ order,
      TRUE ~ class
    ),
    OTU_ID = paste0(ID, "_", OTU_clean),
    OTU_ID = fct_reorder(OTU_ID, percent),
    Fire_Int_Sev = factor(Fire_Int_Sev, levels = fire_levels)
  ) %>%
  ggplot(aes(x = reorder(OTU_ID, percent), y = percent, fill = ID)) +
  geom_col(position = position_dodge()) +
  coord_flip() +
  guides(fill = "none") +
  facet_wrap(~ Fire_Int_Sev,ncol = 1, scales = "free_y") +
  labs(
    x = "",
    y = "",
    tag= '(a)'
  ) +
  theme_classic() +
  theme(
    legend.position = 'top',
    axis.text.x = element_text(hjust = 0.5, size = 25),
    axis.text.y = element_text(hjust = 1, size = 20),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    strip.text = element_text(size = 20),
    plot.tag = element_text(size = 22, hjust = 0.5),
    panel.spacing = unit(1, "lines")
  )

top_freq_soil
