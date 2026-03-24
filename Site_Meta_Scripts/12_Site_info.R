library(tidyverse)


Bag_data<-read.csv('Processed_data/All_Bag_data.csv')

Bag_Sites<-Bag_data%>%
  select(Site,Fire.Severity,Fire.Interval, Vegetation.Class,Total.P,Carbon,Nitrogen,Ammonia_mg_kg,
         Nitrate_mg_kg,Ortho_P_mg_kg,pH,Days_Installed)%>%
  group_by(Site,Fire.Severity,Fire.Interval, Vegetation.Class)%>%
  summarise(round(across(Total.P:Days_Installed, ~ mean(.x, na.rm = TRUE)),digits=2), .groups = "drop")

head(Bag_Sites)
write.csv(Bag_Sites,"Processed_data/12_Sites_Info.csv", row.names=FALSE)


########Extra script below needs to be adjusted########
library(readxl)


Site.info<-read_excel('raw_data/Site.Info.xlsx')


library(sf)
shape <- read_sf("point.dxf", crs= 4326)
#try to convert 

shapes <- shape %>%
  mutate(lon = sf::st_coordinates(.)[,1],
         lat = sf::st_coordinates(.)[,2])

my.sf.point <- st_as_sf(x = shape, 
                        coords = c("longitude", "latitude"),
                        crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Convert to lat long
nc_points <- st_coordinates(st_geometry(st_centroid(shape)))



# Find polygon centroid (This centers the map)
centroid = gCentroid(myshp_proj)
