library(tidyverse)
library(ggplot2)
library(readxl)



#Bray_P####
#Akhter
PO4 <- read_excel("Raw_data/Nutes/Akhter/PO4.xlsx")

Ahkter_P<-PO4%>%
  rename(Site= ABS_plot_id)%>%
  group_by(Site)%>%
  mutate( PO4_mgPL = mean(PO4_mgPL),
    Bray.P.Akhter = PO4_mgPL * 28/4)%>% #mg/kg
select(-sample_id,-absorbance,-survey_phase, -lab_rep)%>%
  distinct()



#Bray_P####
#Sol
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
  filter(Bray.P>=0)%>%
  summarise(Bray.P=mean(Bray.P))


Both_P<-Mean.Bray%>%
  left_join(Ahkter_P%>%
              mutate(Site=as.factor(Site)))
##############
#graph
Both_P%>%
ggplot( aes(x=Site, y=Bray.P))+
  geom_point(aes(alpha=Transect), size= 4)+
  geom_point(aes(x=, y= Bray.P.Akhter), color='red', size = 4)+
  ylab("Bray P (mg/kg)") +
  ylim(0,5)+
  xlab('Site')+
  theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust=.5,size=12),
        axis.title.y = element_text( size= 20) )+
  theme_minimal()



####NH4#####
NH4 <- read_excel("Raw_data/Nutes/Akhter/NH4.xlsx")

NH4.Akhter<-NH4%>%
  rename(Site= ABS_plot_id)%>%
  group_by(Site)%>%
  mutate( NH4_mgNL = mean(NH4_mgNL ),
   NH4.Akhter = NH4_mgNL * 50/5)%>% #mg/kg
select(-Absorbance,-phase,-Time, -lab_rep,-`Manual Dil`,-burnt_status,-Lot)%>%
  distinct()


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
NH4NO3$Result<-NH4NO3$Result * 30/3 #50mL extract for 5g's of soil
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

#removal of extreme outliers
NO3.adj<- NO3[!(NO3$Site == 58 & NO3$Transect == 1), ]
NO3.adj<- NO3.adj[!(NO3.adj$Site == 39 & NO3.adj$Transect == 1), ]


#merge tech replicates
Mean_NH4<-NH4%>%group_by(Site, Transect)%>%
  summarise(NH4=mean(NH4))

#
Both_NH4<-Mean_NH4%>%
  left_join(NH4.Akhter%>%
              mutate(Site=as.factor(Site)))

#graph
Both_NH4%>%
ggplot( aes(x=Site, y= NH4 ))+
  geom_point(aes(y= NH4, alpha=Transect), size= 4)+
  geom_point(aes(y= NH4.Akhter), color='red', size = 4)+
  ylab("NH4 (mg/kg)") +
  xlab('Site')+
  theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust=.5,size=12),
        axis.title.y = element_text( size= 20) )+
  theme_minimal()


#NO3
NO3 <- read_excel("Raw_data/Nutes/Akhter/NO3.xlsx")

NO3.Akhter<-NO3%>%
  rename(Site= ABS_plot_id)%>%
  group_by(Site)%>%
  mutate( NO3_mgNL = mean(NO3_mgNL ),
   NO3.Akhter = NO3_mgNL * 50/5)%>% #mg/kg
select(-NO3_Absorbance,-Time, -lab_rep,-`Manual Dil`,-burnt_status,-Lot, -`Auto Dil`,-NO3_mgNL)%>%
  distinct()

Mean_NO3<-NO3.adj%>%group_by(Site,Transect)%>%
  summarise(NO3=mean(NO3))


Both_NO3<-Mean_NO3%>%
  left_join(NO3.Akhter%>%
              mutate(Site=as.factor(Site)))
#graph
Both_NO3%>%
ggplot( aes(x=Site, y=NO3 ))+
  geom_point(aes(y= NO3, alpha=Transect), size= 4)+
  geom_point(aes(y= NO3.Akhter), color='red', size = 4)+
  ylab("NO3 (mg/kg)") +
  xlab('Site')+
  theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust=.5,size=12),
        axis.title.y = element_text( size= 20) )+
  theme_minimal()


Both_N_and_P<-Both_NH4%>%
  left_join(Both_NO3)%>% select(-NH4_mgNL,-`Auto Dil`)%>%
  left_join(Both_P)%>% select(-PO4_mgPL)

library(writexl)
write_xlsx(Both_N_and_P,'Processed_Data/Ahkter_Sol_N_P_Data.xlsx')
