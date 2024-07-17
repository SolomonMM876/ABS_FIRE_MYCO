


library(tidyverse)
library(ggplot2)
library(readxl)
library(purrr)

Nutrients.Sites.All<-read_excel('Processed_data/Nutrients_Transect_level.xlsx')


myc_data <- read_excel("Raw_data/Bag_data.xlsx", sheet=1)

myc_data$Tube_myc<-as.numeric(myc_data$Tube_myc)
myc_data$myc<-as.numeric(myc_data$myc)
myc_data$beads<-as.factor(myc_data$beads)
myc_data$Tube_ID<-as.factor(myc_data$Tube_ID)


bag_data <- read_excel("Raw_data/Bag_data.xlsx", sheet=2)


bag_data$Location<-as.factor(bag_data$Location)
bag_data$Tube_ID<-as.factor(bag_data$Tube_ID)
bag_data$Site<-as.factor(bag_data$Site)
bag_data$Transect<-as.factor(bag_data$Transect)
bag_data$Nutrient_sub<-as.numeric(bag_data$Nutrient_sub)
bag_data$Bead_weight<-as.numeric(bag_data$Bead_weight)

bag_myc<-left_join(bag_data,myc_data, by='Tube_ID')%>%
  group_by(Site, Transect, Location) %>%
  mutate(Tube_ID = first(Tube_ID),
         #weight of harvested bag
         harvest_w = coalesce(Nutrient_sub, 0) + coalesce(Bead_weight, 0))%>%
  ungroup()

#This is the second round of drying where I subsampled for DNA and CNP measurments
Myc_2nd_round_weight <- read_excel("Raw_data/DW_subsample_DNA_CNP.xlsx", 
                                   sheet = "Sheet2")[,-2] #remove Col called total
Myc_2nd_round_weight$Tube_ID<-as.factor(Myc_2nd_round_weight$Tube_ID)

bag_myc<-left_join(bag_myc, Myc_2nd_round_weight%>%
  mutate(Second_Weight=rowSums(select(., DNA, Chem), na.rm = TRUE))%>%
  select(Tube_ID,Second_Weight), by = 'Tube_ID')


#this is calculated from average bag weight of undeployed bags in bag_weight_vol.R
intial_bag_w<-15.1257*2

#Jeff's math
bag_avg<-bag_myc%>%
  group_by(Site,Transect,Location) %>%
  #remove missing or damaged bags
  filter(is.na(Missing), !Damage %in% c('moderate', 'extreme','major'))%>%
  #remove samples with only 1 rep
  filter(n() == 2) %>%
  mutate(undamaged_harvest_w = sum(harvest_w, na.rm = TRUE))%>%
  group_by(Site,Transect)%>%
  mutate(mean_undam_harvest_w = mean(undamaged_harvest_w))%>%
  select(Site,Transect, mean_undam_harvest_w)%>%
  #join undamaged site avgs
  distinct()%>%
  right_join(bag_myc,by = join_by(Site,Transect))%>%
  #now I repeat what I did above, but for damaged sites as well
  group_by(Site,Transect,Location) %>%
  mutate(harvest_w = sum(harvest_w, na.rm = TRUE))%>%
  #this is my really ugly way of combining the technical replicates
  #I know there has to be a better way, but this just seemed like the best way to avoid data loss
  group_by(Site,Transect,Location)%>%
  mutate(
    mean_undam_harvest_w = first(mean_undam_harvest_w),
    harvest_date = first(harvest_date),
    Nutrient_sub = sum(Nutrient_sub, na.rm = TRUE),
    Bead_weight = sum(Bead_weight, na.rm = TRUE),
    notes.x = paste(notes.x, collapse = " | "),
    roots_in_sample = paste(roots_in_sample, collapse = " | "),
    Damage = paste(Damage, collapse = " | "),
    Missing = paste(Missing, collapse = " | "),
    Before_root.hyphae = first(Before_root.hyphae),
    dirt_in_sample = if_else(any(dirt_in_sample == 'y', na.rm = TRUE), 'y', NA_character_),
    mass_loss = first(mass_loss),
    Tube_w_mg = sum(Tube_w_mg, na.rm = TRUE),
    Tube_myc = sum(Tube_myc, na.rm = TRUE),
    myc = sum(myc, na.rm = TRUE),
    beads = first(beads),
    notes.y = first(notes.y),
    harvest_w = first(harvest_w))%>%
  select(-Rep)%>%
  distinct()%>%
  mutate(initalmass_est = ifelse(grepl('y',Missing) |grepl("moderate|extreme", Damage),
                                 (2*(intial_bag_w * harvest_w) / mean_undam_harvest_w),harvest_w),
         myc_per_bead= (myc / initalmass_est),
         myc_bag_yield_est= myc_per_bead*initalmass_est,
         initalmass_est_avg_all= ((intial_bag_w * harvest_w) / mean_undam_harvest_w),
         myc_per_bead_avg_all= (myc / initalmass_est_avg_all),
         myc_bag_yield_estavg_all= myc_per_bead_avg_all*initalmass_est)%>%
  #now do the same thing to the weights of the myc the second time after weighing
  mutate(Second_Weight_per_bead= (Second_Weight / initalmass_est),
         Second_Weight_bag_yield_est= Second_Weight_per_bead*initalmass_est,
         initalmass_est_avg_all= ((intial_bag_w * harvest_w) / mean_undam_harvest_w),
         Second_Weight_per_bead_avg_all= (Second_Weight / initalmass_est_avg_all),
         Second_Weight_bag_yield_estavg_all= Second_Weight_per_bead_avg_all*initalmass_est)


#CV of bag weights
Site.variation<-bag_avg%>%
  group_by(Site) %>%
  summarise(mean_wt = mean(mean_undam_harvest_w), sd_wt = sd(mean_undam_harvest_w)) %>%
  mutate(cov_wt = sd_wt/mean_wt)%>%
  ungroup()%>%
  mutate(mean_CV=mean(cov_wt))






PROC_VEG_Transect <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                sheet = "Transect.Level_Data")
PROC_VEG_Transect<-PROC_VEG_Transect%>%
  mutate(Site=gsub('ABS00|ABS0', "", Site),
         Transect=gsub('T',"",Transect))

 
Site_Info <- read_excel("Raw_data/Site_Data/Site.Info.xlsx")%>%
  select('1st_Soil_Sample',Bag_Install,Site,Pair)%>%
  mutate(Site=gsub('ABS00|ABS0', "", Site))%>%
           unique()%>%
   group_by(Pair) %>%
  mutate(Site_Pair = paste(sort(unique(Site)), collapse = "-")) %>%
  ungroup()
  

site_data<-list(bag_avg, Nutrients.Sites.All,PROC_VEG_Transect,Site_Info)
#calculate variation within sites in terms of biomass 
Bag_Site<-reduce(site_data,left_join)%>%
  arrange(Transect, Location) %>%
  mutate(Location = factor(Location, levels = sort(unique(Location))))%>%
  group_by(Site)%>%
  mutate(mean_yield = ifelse(myc_bag_yield_est>0, mean(myc_bag_yield_est, na.rm = TRUE), NA),
         sd_yield = sd(myc_bag_yield_est, na.rm = TRUE),
         z_score = (myc_bag_yield_est - mean_yield) / sd_yield,
         is_outlier_Z = abs(z_score) > 2)

rm(bag_data,myc_data,bag_myc,Site_Info,site_data,bag_avg)

#add instillation date
Bag_Site$harvest_date<-as_date(Bag_Site$harvest_date)
Bag_Site$Bag_Install<-as_date(Bag_Site$Bag_Install)


Bag_Site <- Bag_Site %>%
  mutate(Days_Installed = as.integer(harvest_date - Bag_Install))%>%
  arrange(Site,Days_Installed)


mean(Bag_Site$Days_Installed)
sd(Bag_Site$Days_Installed)


Bag_Site$log10_myc_bag_yield_est <- log10(Bag_Site$myc_bag_yield_est + 0.18)
Bag_Site$log10_Second_Weight_bag_yield_est <- log10(Bag_Site$Second_Weight_bag_yield_est + 0.095)
Bag_Site<-Bag_Site%>%
  mutate(Biomass_day=as.numeric(Second_Weight_bag_yield_est/Days_Installed),
         log10_biomass_day= log10(Biomass_day))

hist(Bag_Site$log10_biomass_day)



Bag_Site %>%
  arrange(Second_Weight_bag_yield_est) %>%
  head()%>%
  select(Second_Weight_bag_yield_est)

#hist(log10(.05+Bag_Site$NO3))
Bag_Site$Fire.Interval<-factor(Bag_Site$Fire.Interval, levels = c('Short','Long'))
Bag_Site$Fire.Severity<-factor(Bag_Site$Fire.Severity, levels = c('Low','High'))





library(writexl)
write_xlsx(Bag_Site, path='Processed_data/All_Bag_Site_Info.xlsx')
