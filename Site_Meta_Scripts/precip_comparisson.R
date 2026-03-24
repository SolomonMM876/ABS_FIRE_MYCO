# Load packages
library(raster)
library(sf)
library(chirps)
library(readxl)
library(dplyr)
library(tibble)

# Import site data
Site.info <- read_excel('Raw_data/Site_Data/Site.Info.xlsx')

# Standardize column names
colnames(Site.info)[5] <- 'lat'
colnames(Site.info)[6] <- 'lon'

# Clean site information and assign unique IDs
sites.ID <- Site.info %>%
  distinct(Site, lat, lon, .keep_all = TRUE) %>%
  group_by(Site) %>%
  mutate(ID = cur_group_id())

# Convert to spatial object
sites <- sites.ID %>%
  select(lon, lat) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# ----------------------
# Download CHIRPS data for both periods
# ----------------------

# Period 1: Sep 2023 - Mar 2024
Precip_p1 <- get_chirps(sites, 
                        dates = c("2023-09-01", "2024-03-31"), 
                        server = "ClimateSERV")

# Period 2: Jun 2024 - Nov 2024
Precip_p2 <- get_chirps(sites, 
                        dates = c("2024-06-01", "2024-11-30"), 
                        server = "ClimateSERV")

# Standardize IDs and join with site info
for (df in list(Precip_p1, Precip_p2)) {
  df[[1]] <- as.character(df[[1]])
}
sites.ID$id <- as.factor(sites.ID$ID)

# Calculate mean precipitation per site per period
mean_p1 <- Precip_p1 %>% 
  # group_by(id) %>%
  summarise(mean_precip_p1 = mean(chirps, na.rm = TRUE)) %>% 
  mutate(id=as.factor(id))

mean_p2 <- Precip_p2 %>%
  # group_by(id) %>%
  summarise(mean_precip_p2 = mean(chirps, na.rm = TRUE)) %>% 
  mutate(id=as.factor(id))


# Merge both summaries with site info
precip_comparison <- sites.ID %>%
  select(Site, id, lat, lon) %>%
  left_join(mean_p1, by = "id") %>%
  left_join(mean_p2, by = "id") %>%
  mutate(precip_diff = mean_precip_p2 - mean_precip_p1)

# View result
print(precip_comparison)
