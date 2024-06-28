



library(tidyverse)
library(ggplot2)
library(readxl)



myc_data <- read_excel("~/ABS_FIRE/ABS_FIRE_MYCO/Bag_data.xlsx", sheet=1)

myc_data$Tube_myc<-as.numeric(myc_data$Tube_myc)
myc_data$myc<-as.numeric(myc_data$myc)
myc_data$beads<-as.factor(myc_data$beads)
myc_data$Tube_ID<-as.factor(myc_data$Tube_ID)


bag_data <- read_excel("~/ABS_FIRE/ABS_FIRE_MYCO/Bag_data.xlsx", sheet=2)


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


intial_bag_w<-15.1257*2

#Jeff's math

undamaged_bag_avg<-bag_myc%>%
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
  right_join(bag_myc,by = join_by(Site,Transect))
#now I repeat what I did above, but for damaged sites as well
  group_by(Site,Transect,Location) %>%
  mutate(harvest_w_location = sum(harvest_w, na.rm = TRUE))
  
  
  
###############This is where I am at  
initalmass_est = ((intial_bag_w_location * harvest_w) / undamaged_harvest_w))
  
  
  
  #remove samples that are missing or have moderate|extreme damage
  mutate(undamaged_harvest_w = sum(harvest_w[Missing=='y']))

  
  
  mutate(undamaged_harvest_w = sum(harvest_w[!Missing=='y'|!Damage=="moderate|extreme" ]))
                                             
                                             
                                             
mutate(undamaged_harvest_w = sum(harvest_w[Missing != 'y' & !(Damage %in% c("moderate", "extreme"))]))

,
         




















 
 tmp<-bag_myc%>%
   #what is the average bead weight per site of undamaged bags
   mutate(Combined_weight = coalesce(Nutrient_sub, 0) + coalesce(Bead_weight, 0))%>%
   #remove missing or damaged bags
   filter(is.na(Missing),is.na(Damage)) %>%
   #group by location
   group_by(Site,Transect,Location) %>%
   #remove samples with only 1 rep
   filter(n() == 2) %>%
   mutate(bead_total_weight = sum(Combined_weight, na.rm = TRUE))%>%
   group_by(Site)%>%
   mutate(site_mean_comb_w_bead = mean(bead_total_weight))%>%
   select(Site,Transect site_mean_comb_w_bead)%>%
   #join undamaged site avgs
   distinct()%>%
   right_join(bag_myc,by = join_by(Site))%>%
  group_by(Site,Transect,Location) %>%
  #repeat fnct
   mutate(Combined_weight = coalesce(Nutrient_sub, 0) + coalesce(Bead_weight, 0))%>%
  mutate(bead_total_weight = sum(Combined_weight, na.rm = TRUE))%>%
   #calc ratio in individual bead weight from site average
  mutate( ratio =( site_mean_comb_w_bead/bead_total_weight))%>%
   
#this is my really ugly way of combining the technical replicates
 #I know there has to be a better way, but this just seemed like the best way to avoid data loss
  group_by(Site, Transect, Location) %>%
  mutate(
    site_mean_comb_w_bead = first(site_mean_comb_w_bead),
    harvest_date = first(harvest_date),
    Nutrient_sub = sum(Nutrient_sub, na.rm = TRUE),
    Bead_weight = sum(Bead_weight, na.rm = TRUE),
    Tube_ID = first(Tube_ID),
    notes.x = paste(notes.x, collapse = "|"),
    roots_in_sample = paste(roots_in_sample, collapse = "|"),
    Damage = paste(Damage, collapse = " | "),
    Missing = paste(Missing, collapse = " | "),
    Before_root.hyphae = paste(Before_root.hyphae, collapse = "|"),
    dirt_in_sample = paste(dirt_in_sample, collapse = "|"),
    mass_loss = paste(mass_loss, collapse = "|"),
    Tube_w_mg = sum(Tube_w_mg, na.rm = TRUE),
    Tube_myc = sum(Tube_myc, na.rm = TRUE),
    myc = sum(myc, na.rm = TRUE),
    beads = paste(beads, collapse = "|"),
    notes.y = paste(notes.y, collapse = "|"),
    Combined_weight = sum(Combined_weight, na.rm = TRUE),
    bead_total_weight = first(bead_total_weight),
    ratio = first(ratio)  )%>%
  select(-Rep)%>%
  distinct()%>%
#then I correct for lost myc in locations where a rep was missing or damage was moderate/extreme
   mutate(Corrected_myc = ifelse(grepl('y',Missing) |grepl("moderate|extreme", Damage), myc * ratio, myc))%>%
   mutate(Corrected_myc = ifelse(is.na(Corrected_myc), myc, Corrected_myc))%>%
     select(Damage,Missing,myc,ratio,bead_total_weight,Corrected_myc,site_mean_comb_w_bead)
 

#CV of bag weights
tmp%>%
 group_by(Site) %>%
   summarise(mean_wt = mean(bead_total_weight), sd_wt = sd(bead_total_weight)) %>%
   mutate(cov_wt = sd_wt/mean_wt)%>%
  ungroup()%>%
  mutate(mean_CV=mean(cov_wt))


 #graph combined bag weight 
 tmp%>%
  filter(grepl('y',Missing) |grepl("moderate|extreme", Damage))%>%
  ggplot(aes(x=Site))+
geom_boxplot(aes(color=Transect, y=bead_total_weight))
  


#corrected myc weight
tmp%>%
  ggplot(aes(x=interaction(Site,Transect), y=Corrected_myc))+
           geom_col(aes(fill=Location))

tmp %>%
  ggplot(aes(x = Site, y = Corrected_myc, fill = Location, color=Transect)) +
  geom_col(position = position_dodge(width = 0.9)) +
  labs(x = "Sites", y = "Myc (mg)", fill = "Location") +
  theme_minimal()

bag_myc<-tmp
rm(bag_data,myc_data,tmp)

