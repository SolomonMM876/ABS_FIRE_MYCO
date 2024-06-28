#load libraries and import df
library(raster)
library(dplyr)
library(tibble)
library(readxl)
library(sf)

#import data
Site.info<-read_excel('Raw_data/Updated_Data/ABS.MER.fielddata.Feb.2023_Site.Info_AF.xlsx')

colnames(Site.info)[6]<-'Lat'
colnames(Site.info)[7]<-'Long'

sites <- Site.info[,c('Site', 'Long', 'Lat')] %>%
  distinct(Site, .keep_all = T)%>%
  column_to_rownames(var='Site')


#extract monthly precipitation averages for all locations
#not sure why I cant change file path here, but this is all that works
precip.files <- list.files("Raw_data/world clim/precip", ".tif", full.names=TRUE)
precip <- stack(precip.files)

month <- c("precip_Jan", "precip_Feb", "precip_Mar", "precip_Apr", "precip_May", "precip_Jun", "precip_Jul", "precip_Aug",
           "precip_Sep", "precip_Oct", "precip_Nov", "precip_Dec")
names(precip) <- month

precip.data<-raster::extract(precip,sites)

###################bioclimatic var###
bio.files <- list.files("Raw_data/world clim/bio", ".tif", full.names=TRUE)
bio.var <- stack(bio.files)
#Bio 1 and Bio12 are mean anual temperature and anual precipitation
Temp__Precip_Annual_ <- bio.var[[c(1,12)]]
names(Temp__Precip_Annual_) <- c("Annual_Temp","Annual_Prec")

Temp__Precip_Annual_data<-raster::extract(Temp__Precip_Annual_,sites)


########Tavg######
Tavg.files <- list.files("Raw_data/world clim/tavg", ".tif", full.names=TRUE)
Tavg.var <- stack(Tavg.files)

month <- c("Tavg_Jan", "Tavg_Feb", "Tavg_Mar", "Tavg_Apr", "Tavg_May", "Tavg_Jun", "Tavg_Jul", "Tavg_Aug",
           "Tavg_Sep", "Tavg_Oct", "Tavg_Nov", "Tavg_Dec")
names(Tavg.var) <- month

Tavg_data<-raster::extract(Tavg.var,sites)

###elevation####
elev.files <- list.files("Raw_data/world clim/elevation", ".tif", full.names=TRUE)
elev.var <- stack(elev.files)

elev_data<-raster::extract(elev.var,sites)

###Cbind all extracted data in one df


Site_Env_data<-cbind(sites,precip.data,Temp__Precip_Annual_data, Tavg_data,elev_data)
Site_Env_data<-rownames_to_column(Site_Env_data, var='site')

#remove extra df's
rm(sites,Tavg_data,precip.data,Temp__Precip_Annual_data,elev_data, elav.var,bio.var,Tavg.var, Temp__Precip_Annual_)
