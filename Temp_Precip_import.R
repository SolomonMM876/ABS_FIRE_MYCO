# Install and load necessary packages
library(raster)
library(sf)
library("chirps")
library('readxl')
library('dplyr')
library('tibble')


#change to date
#Wagga_Field_Sites$Date_Burned<-as.Date(as.numeric(Wagga_Field_Sites$Date_Burned),origin = "1899-12-30")


#import data
Site.info<-read_excel('Raw_data/Site_Data/Site.Info.xlsx')


colnames(Site.info)[5]<-'lat'
colnames(Site.info)[6]<-'lon'

sites.ID <- Site.info[,c('Site', 'lon', 'lat')] %>%
  distinct(Site, .keep_all = T)%>%
  group_by(Site)%>%
  mutate(ID = cur_group_id())
  
sites <-as.data.frame(sites.ID[,c('lon','lat')])%>%
  st_as_sf( coords = c("lon", "lat"), crs = 4326) #convert the data frame to an sf object


#extract precipitation data

###this step takes a while
Precip <- get_chirps(sites, dates = c("2023-10-01", "2024-03-10"), server = "ClimateSERV")

#merge df's based on ID
colnames(Precip)[1]<-'ID'
Precip$ID<-as.character(Precip$ID)
sites.ID$ID<-  as.character(sites.ID$ID)

Sites_Precip<-full_join(sites.ID,Precip, by=c('ID','lon','lat'))

#Tmax data from bom Katoomba

Tmax<-read.csv('C:/Users/90957135/OneDrive - Western Sydney University/ABS_FIRE/Tmax-Katoomba/IDCJAC0010_063039_1800_Data.csv')


Tmax$date <- as.Date(with(Tmax, paste(Year, Month, Day,sep="-")), "%Y-%m-%d")

