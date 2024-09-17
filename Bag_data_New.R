# Load required libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(purrr)

# Read in nutrient and site data
Nutrients.Sites.All <- read_excel('Processed_data/Nutrients_Transect_level.xlsx')

# Read in mycorrhizal data from the first sheet of the Excel file
myc_data <- read_excel("Raw_data/Bag_data.xlsx", sheet = 1)
# Convert selected columns to appropriate data types
myc_data <- myc_data %>%
  mutate(
    Tube_myc = as.numeric(Tube_myc),
    myc = as.numeric(myc),
    beads = as.factor(beads),
    Tube_ID = as.factor(Tube_ID)
  )

# Read in bag data from the second sheet of the same Excel file
bag_data <- read_excel("Raw_data/Bag_data.xlsx", sheet = 2)
# Convert selected columns to appropriate data types
bag_data <- bag_data %>%
  mutate(
    Location = as.factor(Location),
    Tube_ID = as.factor(Tube_ID),
    Site = as.factor(Site),
    Transect = as.factor(Transect),
    Nutrient_sub = as.numeric(Nutrient_sub),
    Bead_weight = as.numeric(Bead_weight)
  )

# Consolidate replicate A and B samples at each location into one unit
bag_data <- bag_data %>%
  group_by(Tube_ID) %>%
  mutate(
    harvest_date = first(harvest_date),
    Nutrient_sub = mean(Nutrient_sub, na.rm = TRUE), #so I have this as mean, but for some locations I processed bags separately 
    Bead_weight = sum(Bead_weight, na.rm = TRUE),
    notes = paste(notes, collapse = " | "),
    roots_in_sample = paste(roots_in_sample, collapse = " | "),
    Damage = paste(Damage, collapse = " | "),
    Missing = paste(Missing, collapse = " | "),
    Before_root.hyphae = first(Before_root.hyphae),
    dirt_in_sample = if_else(any(dirt_in_sample == 'y', na.rm = TRUE), 'y', NA_character_),
    mass_loss = first(mass_loss)
  ) %>%
  select(-Rep) %>%
  distinct()
#remove the sample that is a combination of 11-1-16 and 11-2-16
bag_data<-bag_data%>%filter(Tube_ID!= 84)

# Merge bag and mycorrhizal data by Tube_ID
bag_myc <- left_join(bag_data, myc_data, by = 'Tube_ID') %>%
  group_by(Site, Transect, Location) %>%
  mutate(harvest_w = coalesce(Nutrient_sub, 0) + coalesce(Bead_weight, 0)) %>%
  ungroup()

# Read second round weight data for DNA and CNP measurements, excluding "total" column
Myc_2nd_round_weight <- read_excel("Raw_data/DW_subsample_DNA_CNP.xlsx", 
                                   sheet = "Sheet2")[, -2] 
Myc_2nd_round_weight <- Myc_2nd_round_weight %>%
  mutate(Tube_ID = as.factor(Tube_ID))

# Merge second round weight data into bag_myc
bag_myc <- left_join(bag_myc, Myc_2nd_round_weight %>%
                       mutate(Second_Weight = rowSums(select(., DNA, Chem), na.rm = TRUE)) %>%
                       select(Tube_ID, Second_Weight), by = 'Tube_ID')


#this is calculated from average bag weight of undeployed bags in bag_weight_vol.R
initial_bag_w<-15.1257*2
#this has to be the same as Site 29 T1 and added for S29 T2
undamaged_bag_weight_Site_29<- 27.6


  
# Calculate undamaged bag averages
bag_avg <- bag_myc %>%
  #remove missing or damaged bags
  filter(!str_detect(Damage, 'moderate|extreme|major') | str_detect(Missing, 'y')) %>%
  #remove samples with only one bag,locations over 20g because it seems like the weight where the bags arent damaged
  filter(harvest_w > 20) %>%
  group_by(Site, Transect) %>%
  summarise(mean_undam_harvest_w = mean(harvest_w, na.rm = TRUE)) %>%
  ungroup() %>%
  add_row(Site = "29", Transect = "2", mean_undam_harvest_w = undamaged_bag_weight_Site_29) %>%
#now I repeat what I did above, but for damaged sites as well
right_join(bag_myc,by = join_by(Site,Transect))%>%
  group_by(Site,Transect,Location) %>%
  mutate(harvest_w = sum(harvest_w, na.rm = TRUE))%>%
  #correcting the amount of resins found only for dmaged bags
  mutate(initialmass_est = ifelse(harvest_w<20, # I chose 20 because it seems like the weight where the bags arent damaged
                                 (2*(initial_bag_w * harvest_w) / mean_undam_harvest_w),
                                 harvest_w),
        # initialmass_est_all = (mean_undam_harvest_w/harvest_w)*intial_bag_w,
#I am only using the second weight as I believe this to be more accurate without extra sand/beads
         Second_Weight_per_bead= (Second_Weight / harvest_w),
         Second_Weight_est_yield= Second_Weight_per_bead*initialmass_est,
         Second_Weight_est_yield_all_corrected= Second_Weight_per_bead*mean_undam_harvest_w)

#calculations above are as follows:
#calc hypothetical initial mass of beads if bags were damaged
#divide recovered biomass by actual mass of beads that beads were recovered from = mg biomass/g resin
#then multiply above by estimated amount of resins if bags were undamaged




# Load and process vegetation and site information data
PROC_VEG_Transect <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                sheet = "Transect.Level_Data") %>%
  mutate(
    Site = gsub('ABS00|ABS0', "", Site),
    Transect = gsub('T', "", Transect)
  )

Site_Info <- read_excel("Raw_data/Site_Data/Site.Info.xlsx") %>%
  select('1st_Soil_Sample', Bag_Install, Site, Pair) %>%
  mutate(Site = gsub('ABS00|ABS0', "", Site)) %>%
  distinct() %>%
  group_by(Pair) %>%
  mutate(Site_Pair = paste(sort(unique(Site)), collapse = "-")) %>%
  ungroup()

# Merge all site data into one final dataset
site_data <- list(bag_avg, Nutrients.Sites.All, PROC_VEG_Transect, Site_Info)
Bag_Site <- reduce(site_data, left_join)

# Convert dates and calculate days installed
Bag_Site <- Bag_Site %>%
  mutate(
    harvest_date = as_date(harvest_date),
    Bag_Install = as_date(Bag_Install),
    Days_Installed = as.integer(harvest_date - Bag_Install)
  ) %>%
  arrange(Site, Days_Installed)

mean(Bag_Site$Days_Installed)
sd(Bag_Site$Days_Installed)

# Log transformations and biomass calculations
Bag_Site <- Bag_Site %>%
  mutate(
    log10_Second_Weight_bag_yield_est = log10(Second_Weight_est_yield),
    log10_Second_Weight_est_yield_all_corrected = log10(Second_Weight_est_yield_all_corrected),
    Biomass_day = Second_Weight_est_yield / Days_Installed,
    log10_biomass_day = log10(Biomass_day),
    Biomass_day_all_cor = Second_Weight_est_yield_all_corrected / Days_Installed,
    biomass_g_ha_day = Biomass_day_all_cor * (1e+06 / 15),  # Convert to g/ha/day see lines below
    log10_biomass_day_all_cor = log10(Biomass_day_all_cor)
  )
#g/hectare:
#10,000 m^2 ×0.1 m= 1,000 m^3 
#1,000m^3 x (x mg /15cm^3) x (1g/1000mg) x 1000kg/ton x 1000g/kg= g/ha
#(1e+06/15)

# Summary of biomass data
Bag_Site %>%
  summarise(all_mean = mean(biomass_g_ha_day, na.rm = TRUE))

# Plot biomass by site and transect
Bag_Site %>%
  ggplot(aes(x = Site, y = biomass_g_ha_day)) +
  geom_col(aes(fill = Transect), position = "dodge") +
  facet_grid(~Fire.Severity, scales = "free_x")


#Set Factor Order
Bag_Site$Fire.Interval<-factor(Bag_Site$Fire.Interval, levels = c('Short','Long'))
Bag_Site$Fire.Severity<-factor(Bag_Site$Fire.Severity, levels = c('Low','High'))


# CV of bag weights by site
Site.variation <- Bag_Site %>%
  group_by(Site) %>%
  summarise(
    mean_wt = mean(Second_Weight_est_yield_all_corrected, na.rm = TRUE), 
    sd_wt = sd(Second_Weight_est_yield_all_corrected, na.rm = TRUE)
  ) %>%
  mutate(cov_wt = sd_wt / mean_wt) %>%
  ungroup() %>%
  mutate(mean_CV = mean(cov_wt)) %>%
  left_join(PROC_VEG_Transect %>% select(Site, Fire.Interval, Fire.Severity))

# Plot variation by fire interval and biomass
interval_colors <- c("Long" = "darkred", "Short" = "orange")
Site.variation %>%
  ggplot(aes(x = Site, y = cov_wt)) +
  geom_col(aes(fill = Fire.Interval)) +
  facet_grid(~Fire.Interval, scales = "free") +
  labs(x = "Fire Interval", y = "COV of biomass produced") +
  scale_fill_manual(values = interval_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.text.y = element_text(size = 23, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 27, face = "bold"),
    axis.line = element_line(linewidth = 1.5),
    legend.position = 'none'
  )

# Boxplot of Second Weight Yield across sites
Bag_Site %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(Site), y = Second_Weight_est_yield_all_corrected, fill = Transect))


library(writexl)
write_xlsx(Bag_Site, path='Processed_data/All_Bag_Site_Info.xlsx')
