library(readxl)


Site.info<-read_excel('Site.Info.xlsx')


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
