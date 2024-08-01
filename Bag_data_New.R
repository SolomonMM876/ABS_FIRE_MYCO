


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

#this code removes Reps A and B and consolidates the the two bags collected at each location into one unit
bag_data<-bag_data%>%
  group_by(Tube_ID)%>%
  mutate( harvest_date = first(harvest_date),
         Nutrient_sub = mean(Nutrient_sub, na.rm = TRUE),
         Bead_weight = sum(Bead_weight, na.rm = TRUE),
         notes = paste(notes, collapse = " | "),
         roots_in_sample = paste(roots_in_sample, collapse = " | "),
         Damage = paste(Damage, collapse = " | "),
         Missing = paste(Missing, collapse = " | "),
         Before_root.hyphae = first(Before_root.hyphae),
         dirt_in_sample = if_else(any(dirt_in_sample == 'y', na.rm = TRUE), 'y', NA_character_),
         mass_loss = first(mass_loss))%>%
  select(-Rep)%>%
  distinct()
#%>%
 # filter(!Tube_ID==84)




bag_myc<-left_join(bag_data,myc_data, by='Tube_ID')%>%
  group_by(Site, Transect, Location) %>%
  mutate( #weight of harvested bag
    harvest_w = coalesce(Nutrient_sub, 0) + coalesce(Bead_weight, 0))%>%
  ungroup()

#This is the second round of drying where I sub sampled for DNA and CNP measurements
Myc_2nd_round_weight <- read_excel("Raw_data/DW_subsample_DNA_CNP.xlsx", 
                                   sheet = "Sheet2")[,-2] #remove Col called total
Myc_2nd_round_weight$Tube_ID<-as.factor(Myc_2nd_round_weight$Tube_ID)

bag_myc<-left_join(bag_myc, Myc_2nd_round_weight%>%
  mutate(Second_Weight=rowSums(select(., DNA, Chem), na.rm = TRUE))%>%
  select(Tube_ID,Second_Weight), by = 'Tube_ID')


#this is calculated from average bag weight of undeployed bags in bag_weight_vol.R
initial_bag_w<-15.1257*2
#this has to be the same as Site 29 T1 and added for S29 T2
undamaged_bag_weight_Site_29<- 27.6
#Jeff's math
bag_avg<-bag_myc%>%
  group_by(Site,Transect,Location) %>%
  #remove missing or damaged bags
  filter(!str_detect(Damage, 'moderate|extreme|major')|str_detect(Missing,'y'))%>%
  #remove samples with only one bag,locations over 20g because it seems like the weight where the bags arent damaged
  filter(harvest_w>20) %>%
  mutate(undamaged_harvest_w = sum(harvest_w, na.rm = TRUE))%>%
  group_by(Site,Transect)%>%
  mutate(mean_undam_harvest_w = mean(undamaged_harvest_w))%>%
  select(Site,Transect, mean_undam_harvest_w)%>%
  #join undamaged site avgs
  distinct()%>% ungroup()%>%
  #add in undamged bag weight for Site 29 T2 because not enough reps to calculate with above
  add_row(Site = "29", Transect = "2", mean_undam_harvest_w = undamaged_bag_weight_Site_29)%>%
right_join(bag_myc,by = join_by(Site,Transect))%>%
  #now I repeat what I did above, but for damaged sites as well
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
Bag_Site<-reduce(site_data,left_join)

rm(bag_data,myc_data,bag_myc,Site_Info,site_data,bag_avg)

#add instillation date
Bag_Site$harvest_date<-as_date(Bag_Site$harvest_date)
Bag_Site$Bag_Install<-as_date(Bag_Site$Bag_Install)


Bag_Site <- Bag_Site %>%
  mutate(Days_Installed = as.integer(harvest_date - Bag_Install))%>%
  arrange(Site,Days_Installed)


mean(Bag_Site$Days_Installed)
sd(Bag_Site$Days_Installed)


Bag_Site$log10_Second_Weight_bag_yield_est <- log10(Bag_Site$Second_Weight_est_yield)
Bag_Site$log10_Second_Weight_est_yield_all_corrected <- log10(Bag_Site$Second_Weight_est_yield_all_corrected)

#g/hectare:
#10,000 m^2 ×0.1 m= 1,000 m^3 
#1,000m^3 x (x mg /15cm^3) x (1g/1000mg) x 1000kg/ton x 1000g/kg= g/ha
(1e+06/15)


Bag_Site<-Bag_Site%>%
  mutate(Biomass_day=as.numeric(Second_Weight_est_yield/Days_Installed),
         log10_biomass_day= log10(Biomass_day),
         Biomass_day_all_cor=as.numeric(Second_Weight_est_yield_all_corrected/Days_Installed),
         biomass_g_ha_day= Biomass_day_all_cor*(1e+06/15),#convert to g/ha/day 
         log10_biomass_day_all_cor= log10(Biomass_day_all_cor))

Bag_Site%>%
  summarise(all_mean= mean(biomass_g_ha_day, na.rm=TRUE))
         


Bag_Site%>%
  ggplot(aes(x=Site,y=biomass_g_ha_day))+
           geom_col(aes(fill=Transect),position = "dodge")+
  facet_grid(~Fire.Severity,scales = "free_x")



Bag_Site %>%
  arrange(Second_Weight_est_yield) %>%
  head()%>%
  select(Second_Weight_est_yield)

#hist(log10(.05+Bag_Site$NO3))
Bag_Site$Fire.Interval<-factor(Bag_Site$Fire.Interval, levels = c('Short','Long'))
Bag_Site$Fire.Severity<-factor(Bag_Site$Fire.Severity, levels = c('Low','High'))





library(writexl)
write_xlsx(Bag_Site, path='Processed_data/All_Bag_Site_Info.xlsx')
