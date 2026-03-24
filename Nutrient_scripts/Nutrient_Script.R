####Soil Nutrient Data. Decomp####

library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)
#library(splitstackshape)
library(readxl)
library(car)
library(ggrepel)




###Site Descriptions####
Site_Info<-read_excel("Raw_data/Sequence/fastfa_Site.Info.xlsx")
colnames(Site_Info)[3]<-'Site'
colnames(Site_Info)[4]<-'Transect'
Site_Info$Transect<-as.factor(Site_Info$Transect)
Site_Info$Site= sub(c('ABS00'),'',Site_Info$Site)
Site_Info$Site= sub(c('ABS0'), '',Site_Info$Site)

#Bray_P####
#Load data
Bray_P <- lapply(c('Raw_Data/Nutes/23-06-07 solomon bray -P samples lot1 2023.csv',
                   'Raw_Data/Nutes/23-06-08 solomon bray-P samples lot2.csv',
                    'Raw_Data/Nutes/23-06-09 solomon bray-p samples lot3.csv'),read.csv)
#Clean
Bray_P<-bind_rows(Bray_P, .id = "sample.ID")
Bray_P$Sample.ID<-sub(c("b"),"",Bray_P$Sample.ID)
Bray_P$Sample.ID<-sub(c("a"),"",Bray_P$Sample.ID)
Bray_P$Sample.ID<-sub(c(" "),"-",Bray_P$Sample.ID)
colnames(Bray_P)[2]<-'Transect'
colnames(Bray_P)[5]<-'Bray.P'
Bray_P$Bray.P<-Bray_P$Bray.P * 50/3.57 #50mL extract for 3.57g's of soil
Bray_P$Units<-'mg/kg'
Bray_P<-Bray_P %>% separate(Transect,c('Site','Transect'),extra = "merge")
Bray_P$Site<-as.factor(as.numeric(Bray_P$Site))
Bray_P$Bray.P<-as.numeric(Bray_P$Bray.P)
Bray_P<-Bray_P %>%drop_na(Site)


#merge tech replicates
Mean.Bray<-Bray_P%>%group_by(Site,Transect)%>%
  mutate(Bray.P = ifelse(Bray.P<0, 0.05, Bray.P))%>%
  summarise(Bray.P=mean(Bray.P))

Bray_P.Site<-left_join(Mean.Bray,Site_Info, by = join_by(Site,Transect))
Bray_All<-left_join(Bray_P,Site_Info, by = join_by(Site,Transect))

#Anova
ggplot(Bray_P.Site, aes( x=Bray.P))+
  geom_histogram()

qqPlot(Bray_P.Site$Bray.P)

Bray_Aov<-lm(Bray.P~Interval+Severity, data=Bray_P.Site)
summary(Bray_Aov)

All.Bray.P<-left_join(Site_Info,Mean.Bray, by = join_by(Site,Transect))



##############
#graph
Bray_All%>%
  mutate(Site = fct_reorder(Site,Bray.P))%>%
ggplot( aes(x=Site, y=Bray.P, alpha=Transect))+
  geom_point()+
  facet_wrap(~Interval,scales = "free_x")+
  ylab("Bray P (mg/kg)") +
  xlab('Site (Transect 1+2)')+
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust=.5))


#setwd('graphs')
#ggsave(filename=paste('Bray P_Interval.png'), width = 25, height = 15, dpi = 800, units = "cm", device='png')




















#N####
#setwd("\\\\ad.uws.edu.au/dfshare/HomesCMB$/90957135/My Documents/ABS_FIRE/Fire_Decomp_Nuts")
#import
NH4NO3<-lapply(c('Raw_Data/Nutes/23-06-13 solomon NH4NO3 LOT1.csv',
                 'Raw_Data/Nutes/23-06-13 solomon NH4NO3 LOT2.csv',
                 'Raw_Data/Nutes/23-06-15 solomon nh4no3 lot3 .csv',
                 'Raw_Data/Nutes/23-06-16 solomon NH4NO3 LOT4 LAST.csv'),read.csv)

#clean and organize
NH4NO3<-bind_rows(NH4NO3, .id = "sample.ID")
NH4NO3$Sample.ID<-sub(c("a"),"",NH4NO3$Sample.ID)
NH4NO3$Sample.ID<-sub(c("b"),"",NH4NO3$Sample.ID)
NH4NO3$Sample.ID<-sub(c(" "),"-",NH4NO3$Sample.ID)
#convert units
NH4NO3$Result<-NH4NO3$Result * 50/5 #50mL extract for 5g's of soil
NH4NO3$Units<-'mg/kg'
colnames(NH4NO3)[2]<-'Transect'
NH4NO3<-NH4NO3 %>% separate(Transect,c('Site','Transect'),extra = "merge")
NH4NO3$Site<-as.factor(as.numeric(NH4NO3$Site))
NH4NO3<-NH4NO3 %>%drop_na(Site)
str(NH4NO3)
#remove sites that are not part of 30 sites
#54,55,61,63

New.NH4NO3<-subset(NH4NO3,!Site %in% c(54,55,61,63))



#split to two dfs
NH4<-New.NH4NO3[New.NH4NO3$Test.Name =="Ammonia 2.0",]
colnames(NH4)[6]<-'NH4'

NO3<-New.NH4NO3[New.NH4NO3$Test.Name =="Nitrate 2",]
colnames(NO3)[6]<-'NO3'

#merge tech replicates
Mean_NH4<-NH4%>%group_by(Site, Transect)%>%
  summarise(NH4=mean(NH4))


All.Mean.NH4<-left_join(Site_Info,Mean_NH4, by = join_by(Site,Transect))


NH4.Site<-left_join(Mean_NH4,Site_Info, by = join_by(Site,Transect))
NH4_All<-left_join(NH4,Site_Info, by = join_by(Site,Transect))

#merge tech replicates
Mean_NO3<-NO3%>%group_by(Site,Transect)%>%
  summarise(NO3=mean(NO3))


All.Mean.NO3<-left_join(Site_Info,Mean_NO3, by = join_by(Site,Transect))



NO3.Site<-left_join(Mean_NO3,Site_Info, by = join_by(Site,Transect))
NO3_All<-left_join(NO3,Site_Info, by = join_by(Site,Transect))

#Anova
ggplot(NH4.Site, aes( x=NH4))+
  geom_histogram()

qqPlot(NH4.Site$NH4)

NH4_Aov<-aov(NH4~Interval+Severity, data=NH4.Site)
summary(NH4_Aov)

##############
#graph
#NH4
NH4_All%>%
  mutate(Site = fct_reorder(Site,NH4))%>%
ggplot( aes(x=Site, y=NH4, alpha=Transect))+
  geom_point()+
  facet_wrap(~Severity,scales = "free_x")+
  ylab("NH4 (mg /kg)") +
  xlab('Site (Transect 1+2)')+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0))

#setwd('graphs')
#ggsave(filename=paste('NH4_Severity.png'), width = 40, height = 15, dpi = 800, units = "cm", device='png')

#NO3



#Anova
ggplot(NO3.Site, aes( x=NO3))+
  geom_histogram()

qqPlot(NO3.Site$NO3)

NO3_Aov<-aov(NO3~Interval+Severity, data=NO3.Site)
summary(NO3_Aov)

NO3_All%>%
  mutate(Site = fct_reorder(Site,NO3))%>%
ggplot( aes(x=Site, y=NO3, alpha=Transect))+
  geom_point()+
  facet_wrap(~Interval,scales = "free_x")+
  ylab("NO3 (mg /kg)")+
  xlab('Site (Transect 1+2)')+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0))

#setwd('graphs')
#ggsave(filename=paste('NO3_Interval.png'), width = 40, height = 15, dpi = 800, units = "cm", device='png')




#Total_P####
#setwd("\\\\ad.uws.edu.au/dfshare/HomesCMB$/90957135/My Documents/ABS_FIRE/Fire_Decomp_Nuts")
#Load Data
Total_P<-read.csv("Raw_Data/Nutes/CSBP_Total_P.csv")
names(Total_P)[6]<-'Transect'

#clean and organize
colnames(Total_P)<-Total_P[3,]
ncol(Total_P)
colnames(Total_P)[8]<-"Total.P"
colnames(Total_P)[6]<-"Transect"
Total_P<-Total_P[-c(1,2,3,65),]
Total_P$Total.P<-as.numeric(Total_P$Total.P)
Total_P$Transect<-as.character(Total_P$Transect)

#remove empty columns 
Total_P<-Total_P[-c(2,4)]

Total_P<-separate(Total_P,Transect,into=c('Site','Transect'),sep='-')
Total_P$Site<-as.factor(as.numeric(Total_P$Site))


Mean.T.P<-Total_P%>%group_by(Site,Transect)%>%
  mutate(Total.P=mean(Total.P))%>%
  select(Total.P,Site,Transect)

Mean.T.P<-Mean.T.P%>%
  group_by(Site,Transect)%>%
  summarise(Total.P=mean(Total.P))


All.Mean.Total.P<-left_join(Site_Info,Mean.T.P, by = join_by(Site,Transect))



Total_P.Site<-left_join(Mean.T.P,Site_Info, by = join_by(Site,Transect))
Total_P_All<-left_join(Total_P,Site_Info, by = join_by(Site,Transect))

#Anova
ggplot(Total_P.Site, aes( x=Total.P))+
  geom_histogram()

qqPlot(Total_P.Site$Total.P)

Total_P_Aov<-aov(Total.P~Interval+Severity, data=Total_P.Site)
summary(Total_P_Aov)

##############
#graph
Total_P_All%>%
  mutate(Site = fct_reorder(Site,Total.P))%>%
  ggplot( aes(x=Site, y=Total.P, alpha=Transect))+
  geom_point()+
  facet_wrap(~Severity,scales = "free_x")+
  ylab("Total P (mg/kg)") +
  xlab('Site (Transect 1+2)')+
  theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust=.5))

#setwd('graphs')
#ggsave(filename=paste('Total_P_Severity.png'), width = 25, height = 15, dpi = 800, units = "cm", device='png')









#Total_C+N####
#setwd("\\\\ad.uws.edu.au/dfshare/HomesCMB$/90957135/My Documents/ABS_FIRE/Fire_Decomp_Nuts")
#Load Data
Total_CN<-read.csv("Raw_Data/Nutes/CN SOIL SAMPLES_2.0.csv")

colnames(Total_CN)[1]<-'Transect'
colnames(Total_CN)[5]<-'Carbon'
colnames(Total_CN)[6]<-'Nitrogen'
Total_CN<-Total_CN%>% separate(Transect, c('Site','Transect'),extra = "merge")
Total_CN$Site<-as.factor(as.numeric(Total_CN$Site))
Total_CN$Transect<-sub(c("T"),"",Total_CN$Transect)
Total_CN$Carbon<-as.numeric(Total_CN$Carbon)
Total_CN$Nitrogen<-as.numeric(Total_CN$Nitrogen)



Mean.T.C.N<-Total_CN%>%
  mutate(Transect = gsub('A|B','',Transect))%>%
  filter(!is.na(Transect))%>%
  group_by(Site,Transect)%>%
  summarise(Carbon=mean(Carbon),
            Nitrogen=mean(Nitrogen))

Total_CN.Site<-left_join(Mean.T.C.N,Site_Info, by = join_by(Site,Transect))
Total_CN_All<-left_join(Total_CN,Site_Info, by = join_by(Site,Transect))



















#XRF P and S####
#setwd("\\\\ad.uws.edu.au/dfshare/HomesCMB$/90957135/My Documents/ABS_FIRE/Fire_Decomp_Nuts")
#load data
XRF_S_P<-read_tsv('Raw_Data/Nutes/All Sites.edit.txt', show_col_types = FALSE)
colnames(XRF_S_P)[2]<-'Transect'
colnames(XRF_S_P)[6]<-'XRF_P'
colnames(XRF_S_P)[8]<-'XRF_S'
XRF_S_P<-XRF_S_P %>% separate(Transect,c('Site','Transect'))
#remove SD and Avg
XRF_S_P<-XRF_S_P[-c(65,66),]
XRF_S_P$Site<-as.factor(as.numeric(XRF_S_P$Site))


#Labels for x axis
#Makes it easier to read x axis
XRF_S_P.Lab<-XRF_S_P%>%distinct(Transect,Site)
Site.XRF_S_P<-XRF_S_P.Lab[order(XRF_S_P.Lab$Site),]$Site

ggplot(XRF_S_P, aes(x=interaction(Transect,Site), y=XRF_P, alpha=Transect))+
  geom_point()+
  scale_x_discrete(label=Site.XRF_S_P)+
  ylab("P XRF (mg /kg)")+
  xlab('Site (Transect 1+2)')+
  theme(axis.text.x = element_text(angle = -45, vjust = -.4, hjust=1))
#setwd('graphs')
#ggsave(filename=paste('XRF_P.png'), width = 25, height = 15, dpi = 800, units = "cm", device='png')

ggplot(XRF_S_P, aes(x=interaction(Transect,Site), y=XRF_S, alpha=Transect))+
  geom_point()+
  scale_x_discrete(label=Site.XRF_S_P)+
  ylab("S XRF (mg /kg)")+
  xlab('Site (Transect 1+2)')+
  theme(axis.text.x = element_text(angle = -45, vjust = -.4, hjust=.5))
#setwd('graphs')
#ggsave(filename=paste('XRF_S.png'), width = 25, height = 15, dpi = 800, units = "cm", device='png')



###Compare XRF to Total P
Ps<-left_join(XRF_S_P,Total_P, by = join_by(Site,Transect), unmatched = "drop",relationship = "many-to-many")



lm_model <- lm(XRF_P~ Total.P, data = Ps)


# Extract R-squared value
r_squared <- summary(lm_model)$r.squared

library(ggpmisc)
Ps%>%
  arrange(XRF_P)%>%
ggplot((aes(x=Total.P, y=XRF_P))) +
  stat_poly_line() +
  stat_poly_eq() +
  geom_point()
#setwd('graphs')
#ggsave(filename=paste('XRF P vs Lab P.png'), width = 40, height = 15, dpi = 800, units = "cm", device='png')

##############


######Nutrient data#######




df_list_All<-list(Mean.Bray, Mean_NH4, Mean_NO3,Mean.T.P, Mean.T.C.N ) 

Nutrients.Sites.All<-Reduce(function(x, y) merge(x, y, all=T), df_list_All)


Nutrients.Site<-Nutrients.Sites.All%>%
  select(-Transect)%>%
  group_by(Site)%>%
  mutate(across(c(Bray.P:Nitrogen),mean, na.rm= TRUE))%>%
  unique()


# List all objects in the current environment
all_objects <- ls()

# Specify objects to keep
objects_to_keep <- c("Nutrients.Sites.All", "Site_Info","Nutrients.Site")

# Identify objects to remove (all objects except those to keep)
objects_to_remove <- setdiff(all_objects, objects_to_keep)

# Remove identified objects
rm(list = objects_to_remove)


library(writexl)
write_xlsx(Nutrients.Site, path='Processed_data/Nutrients_Site_level.xlsx')
write_xlsx(Nutrients.Sites.All, path='Processed_data/Nutrients_Transect_level.xlsx')


