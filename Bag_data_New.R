# Load required libraries
library(tidyverse)
library(ggplot2)
library(readxl)

# Read in mycorrhizal data from the first sheet of the Excel file
myc_data <- read_excel("Raw_data/Bag_data.xlsx", sheet = 1)
# Convert selected columns to appropriate data types
myc_data <- myc_data %>%
  mutate(
    Tube_myc_w = as.numeric(Tube_myc),
    myc_w = as.numeric(myc),
    resins = as.factor(beads),
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
    Res_nute_sub_w = as.numeric(Nutrient_sub),
    Res_remain_w = as.numeric(Bead_weight)
  )

# Consolidate replicate A and B samples at each location into one unit
bag_data <- bag_data %>%
  group_by(Tube_ID,Site,Transect,Location) %>%
  summarise(
    harvest_date = first(harvest_date),
    Res_nute_sub_w = sum(Res_nute_sub_w, na.rm = TRUE),#this is a sum here now
    Res_remain_w = sum(Res_remain_w, na.rm = TRUE),
    notes = paste(notes, collapse = " | "),
    roots_in_sample = paste(roots_in_sample, collapse = " | "),
    Damage = paste(Damage, collapse = " | "),
    Missing = paste(Missing, collapse = " | "),
    Before_root.hyphae = first(Before_root.hyphae),
    dirt_in_sample = if_else(any(dirt_in_sample == 'y', na.rm = TRUE), 'y', NA_character_),
    mass_loss = first(mass_loss)
  ) %>%
  #remove the sample that is a combination of 11-1-16 and 11-2-16 because of processing error
  filter(Tube_ID!= 84)
#Samples 29-1-33 and samples 56-1-16 had both bags damaged so there was no biomass collected from these samples

library(writexl)
write.csv(bag_data, file='Processed_data/bag_data.csv')


# Merge bag and mycorrhizal data by Tube_ID
bag_myc <- left_join(bag_data, myc_data, by = 'Tube_ID') %>%
  group_by(Site, Transect, Location) %>%
  mutate(Res_total_bag= Res_nute_sub_w+Res_remain_w)%>%
  ungroup()

# Read second round weight data for DNA and CNP measurements, excluding "total" column
Myc_2nd_round_weight <- read_excel("Raw_data/Stoich/DW_subsample_DNA_CNP.xlsx", 
                                   sheet = "Sheet2")[, -2]  %>%
  mutate(Tube_ID = as.factor(Tube_ID))

# Merge second round weight data into bag_myc
bag_myc <- left_join(bag_myc, Myc_2nd_round_weight %>%
                       mutate(myc_2nd_w = rowSums(select(., DNA, Chem), na.rm = TRUE)) %>%
                       select(Tube_ID, myc_2nd_w), by = 'Tube_ID')

dat_temp<-bag_myc%>%group_by(Site,Transect,Location) %>%
  summarise(Res_total_Location = sum(Res_total_bag, na.rm = TRUE),
            myc_2nd_w=sum(myc_2nd_w),
            Res_remain_w=sum(Res_remain_w))


#this is calculated from average bag weight of undeployed bags in bag_weight_vol.R
initial_bag_w<-15.1257*2
#this has to be the same as Site 29 T1 and added for S29 T2
undamaged_bag_weight_Site_29<- 27.6
  
# Calculate undamaged bag averages
bag_avg_undamaged <- bag_myc %>%
  #remove missing or damaged bags
  filter(!str_detect(Damage, 'moderate|extreme|major') | str_detect(Missing, 'y')) %>%
  #remove samples with only one bag,locations over 20g because it seems like the weight where the bags arent damaged
  filter(Res_total_bag > 20) %>%
  group_by(Site, Transect) %>%
  summarise(mean_undam_resin_w = mean(Res_total_bag, na.rm = TRUE)) %>%
  ungroup() %>%
  add_row(Site = "29", Transect = "2", mean_undam_resin_w = undamaged_bag_weight_Site_29)

#here I am calculating the amount of biomass of myc from each location accounting for variation in soil moisture from each site
corrected_myc<-left_join(dat_temp,bag_avg_undamaged)%>%
  #(initial bag weight saturated * amount Resin collected per Location)/ avg weight of resins per location at each transect
  mutate(resin_mass_est = (initial_bag_w * Res_total_Location) / mean_undam_resin_w,
         myc_2nd_w_per_bead= (myc_2nd_w / resin_mass_est),
         myc_2nd_w_est_yield= myc_2nd_w_per_bead*initial_bag_w)


#Site info data
Site_Info <- read_excel("Raw_data/Site_Data/Site.Info.xlsx") %>%
  select('1st_Soil_Sample', Bag_Install, Site, Pair, Latitude,Longitude) %>%
  mutate(Site = gsub('ABS00|ABS0', "", Site)) %>%
  distinct() %>%
  group_by(Pair) %>%
  mutate(Site_Pair = paste(sort(unique(Site)), collapse = "-")) %>%
  ungroup()


# Merge all site data into one final dataset
site_data <- left_join(bag_data%>%select(Site,Transect,Location,harvest_date), Site_Info)%>%
# Convert dates and calculate days installed
  mutate(
    harvest_date = as_date(harvest_date),
    Bag_Install = as_date(Bag_Install),
    Days_Installed = as.integer(harvest_date - Bag_Install)
  ) %>%
  arrange(Site, Days_Installed)

mean(site_data$Days_Installed)
sd(site_data$Days_Installed)

# Log transformations and biomass calculations
Bag_Site <- left_join(corrected_myc,site_data) %>%
  mutate(
    log10_myc_2nd_w_est_yield = log10(myc_2nd_w_est_yield),
    Biomass_day = myc_2nd_w_est_yield / Days_Installed,
    log10_biomass_day = log10(Biomass_day),
    biomass_g_ha_day = Biomass_day * (1e+06 / 15),  # Convert to g/ha/day see lines below
  )
#g/hectare:
#10,000 m^2 ×0.1 m= 1,000 m^3 
#1,000m^3 x (x mg /15cm^3) x (1g/1000mg) x 1000kg/ton x 1000g/kg= g/ha
#(1e+06/15)

# Summary of biomass data
Bag_Site %>%
  summarise(all_mean = mean(biomass_g_ha_day, na.rm = TRUE))


# Load and process vegetation and site information data
VEG_Transect <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                sheet = "Transect.Level_Data") %>%
  mutate(
    Site = gsub('ABS00|ABS0', "", Site),
    Transect = gsub('T', "", Transect)
  )
# Read in nutrient data
Nutrients.Sites.T <- read_excel('Processed_data/Nutrients_Transect_level.xlsx')

Bag_Site<-Bag_Site%>%
  left_join(VEG_Transect)%>%
  left_join(Nutrients.Sites.T)

#Set Factor Order
Bag_Site$Fire.Interval<-factor(Bag_Site$Fire.Interval, levels = c('Short','Long'))
Bag_Site$Fire.Severity<-factor(Bag_Site$Fire.Severity, levels = c('Low','High'))

# Plot biomass by site and transect
Bag_Site %>%
  ggplot(aes(x = Site, y = biomass_g_ha_day)) +
  geom_col(aes(fill = Transect), position = "dodge") +
  facet_grid(~Fire.Interval, scales = "free_x")



# CV of bag weights by site
Site.variation <- Bag_Site %>%
  group_by(Site) %>%
  summarise(
    mean_wt = mean(myc_2nd_w_est_yield, na.rm = TRUE), 
    sd_wt = sd(myc_2nd_w_est_yield, na.rm = TRUE)
  ) %>%
  mutate(cov_wt = sd_wt / mean_wt) %>%
  ungroup() %>%
  mutate(mean_CV = mean(cov_wt)) %>%
  left_join(VEG_Transect %>% select(Site, Fire.Interval, Fire.Severity))

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
  geom_boxplot(aes(x = as.factor(Site), y = myc_2nd_w_est_yield, fill = Transect))



write.csv(Bag_Site, file='Processed_data/All_Bag_Site_Info.csv',row.names = FALSE)

