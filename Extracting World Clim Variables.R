
#Extracting World Clim variables

# can get precipition and temperature by querying worldclim within R
# see https://www.worldclim.org/data/worldclim21.html for more info
#load once then reload from saved file
#prec <- getData('worldclim', var='prec', res=2.5)
#writeRaster(prec, "Raw_data/clim/prec.tif", format = "GTiff", overwrite = TRUE)
prec <- raster("Raw_data/clim/prec.tif")

#tavg <- getData('worldclim', var='tmean', res=2.5)
#writeRaster(tavg, "Raw_data/clim/tavg.tif", format = "GTiff", overwrite = TRUE)
tavg <- raster("Raw_data/clim/tavg.tif")

#elav<-raster::getData('alt',country='AUS')
#writeRaster(elav, "Raw_data/clim/elav.tif", format = "GTiff", overwrite = TRUE)
tavg <- raster("Raw_data/clim/elav.tif")


# extract values at sampled coordinates and add to metadata table
meta$prec <- apply(extract(prec, meta[, c('Longitude', 'Latitude')]), 1, sum)
meta$tmean <- apply(extract(tavg, meta[, c('Longitude', 'Latitude')]), 1, mean)
meta$tmean <- meta$tmean/10  
