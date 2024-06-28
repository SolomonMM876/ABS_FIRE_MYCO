
library(dplyr)
library(vegan)

Bag_Site<-read_excel('~/ABS_FIRE/ABS_FIRE_MYCO/Processed_data/All_Bag_Site_Info.xlsx')
Nutrients.Sites.All<-read_excel('~/ABS_FIRE/ABS_FIRE_MYCO/Processed_data/Nutrients_Site_level.xlsx')
PROC_VEG_Site <- read_excel("~/ABS_FIRE/ABS_FIRE_MYCO/Raw_data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                sheet = "Site.Level_Data")%>%
      mutate(Site=gsub('ABS00|ABS0', "", Site))

names(PROC_VEG_Site)[6]<-'Second_Fire_yr'  
names(PROC_VEG_Site)[7]<-'Third_Fire_yr'  
PROC_VEG_Site$Third_Fire_yr<-as.numeric(PROC_VEG_Site$Third_Fire_yr)



#drop rows with 0 biomass
x <- which(rowSums(select(Bag_Site, myc)) == 0)
Bag_Site[x, ]
Bag_Site <- Bag_Site[-x, ]  # drop those rows

#Mycoriz data
Myc_data<-Bag_Site%>%
  select(Site,Transect,Location,log10_myc_bag_yield_est,Days_Installed,Pair)


#Nutrient df rows=Site
Site.Nutrients<-data.frame(Nutrients.Sites.All, row.names=1)
#NH4 and NO3 highly correlated and need to pick one

#Veg data df row= Site
#remove extra sites that are not part of the analysis
#Site.Veg.Factors has site info as well
Site.Veg.Factors <- PROC_VEG_Site %>%
  semi_join(Nutrients.Sites.All, by = "Site")%>%
  mutate(First_Interval = as.integer(Most.Recent.Fire_Year - Second_Fire_yr),
         Second_Interval = if_else(is.na(Third_Fire_yr)==FALSE, as.integer(Second_Fire_yr - Third_Fire_yr),NA),
         Interval_Avg= (First_Interval+Second_Interval)/2)%>%
  select(-c(Vegetation.Class,Most.Recent.Fire_Year,Second_Fire_yr,Third_Fire_yr,First_Interval,Second_Interval,
            Sample.Date_Fieldwork,Fire.Type))%>%
  data.frame(row.names = 1)

Site.Veg<-Site.Veg.Factors%>%
  select(Tree.Basal.Area_m2:Interval_Avg)







#Nute.pca <- rda( Nutrient.num_r~ Severity + Interval, data=Treat.w.nutes, scale=TRUE)
Nute.pca <- rda( Nutrient.num_r, data=Treat.w.nutes, scale=TRUE)

#Nute.pca<-rda(Nutrient.num, scale=TRUE)
plot(Nute.pca)    
summary(Nute.pca)


Nute.spe2.sc <- scores(Nute.pca, tidy=T)
Nute.spe2.sc.ar <- Nute.spe2.sc %>% 
  filter(score=='species')
Nute.spe2.sc.ar <- tibble::rownames_to_column(Nute.spe2.sc.ar)


Nute_Sites <- Nute.spe2.sc %>%
  filter(score=='sites') %>% 
  rename(Site=label) %>% 
  left_join(Site_Info)

p<-ggplot(Nute_Sites, aes(x=PC1, y=PC2, colour=Severity, shape=Interval, label=Site )) + 
  geom_point(size=5)+ 
  geom_segment(data=Nute.spe2.sc.ar,inherit.aes = FALSE,
               aes(x=0,y=0, xend=PC1, yend=PC2, group=rowname),
               arrow = arrow(type = "closed",length=unit(4,'mm')),
               color= 'black') +
  geom_text(colour='black',size=4)+
  geom_text_repel(data= Nute.spe2.sc.ar, inherit.aes = FALSE,
                  aes(x= PC1,y=PC2,label=rowname))+
  labs(x = "PC1 55.8%", y = "PC2 16.6%")+ theme_minimal()

p
#setwd('graphs')
#ggsave(filename=paste('Nute_PCA_num.png'), width = 25, height = 15, dpi = 800, units = "cm", device='png')


plotly::ggplotly(p)

###Veg data####
#Load Data

Veg_Cover<-read_excel("Raw_data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", sheet="Site.Level_Data")
head(Veg_Cover)
colnames(Veg_Cover)[1]<-'Site'
Veg_Cover$Site= sub(c('ABS00'),'',Veg_Cover$Site)
Veg_Cover$Site= sub(c('ABS0'), '',Veg_Cover$Site)
head(Veg_Cover)
colnames(Veg_Cover)[4]<-'Interval'
colnames(Veg_Cover)[8]<-'Severity'
colnames(Veg_Cover)[6]<-'Second_Most_Receant_Fire_Year'
colnames(Veg_Cover)[7]<-'Third_Most_Receant_Fire_Year'

#now join both df so only decomp sites are present in df

Veg_Cover_Decomp<-left_join(Site_Info[3],Veg_Cover)
Veg_Cover<-as.data.frame(Veg_Cover_Decomp)
#remove col with characters
Veg_Cover.num <- Veg_Cover_Decomp %>%
  subset(select=-c( Vegetation.Class:Sample.Date_Fieldwork))


colnames(Veg_Cover.num)
Veg_Cover.num <- mutate_all(Veg_Cover.num, function(x) as.numeric(as.character(x)))
Veg_Cover.num<-as.data.frame(Veg_Cover.num)






Veg_Cover_r<-data.frame(Veg_Cover, row.names=1)
Veg_Cover_num.r<-data.frame(Veg_Cover.num, row.names=1)




Veg_Nute.pca <- rda( Veg_Cover_num.r~ Severity + Interval, data=Veg_Cover_r, scale=TRUE)
Veg.pca <- rda( Veg_Cover_num.r, data=Veg_Cover_r, scale=TRUE)
###Making PCA#

plot(Veg.pca)    
summary(Veg.pca)


Veg.spe2.sc <- scores(Veg.pca, tidy=T)
Veg.spe2.sc.ar <- Veg.spe2.sc %>% 
  filter(score=='species')
Veg.spe2.sc.ar <- tibble::rownames_to_column(Veg.spe2.sc.ar)


Veg_Sites <- Veg.spe2.sc %>%
  filter(score=='sites') %>% 
  rename(Site=label) %>% 
  left_join(Site_Info)

#Veg_Sites<-na.omit(Veg_Sites)


p<-ggplot(Veg_Sites, aes(x=PC1, y=PC2, colour=Severity, shape=Interval, label=Site )) + 
  geom_point(size=5)+ 
  geom_segment(data=Veg.spe2.sc.ar,inherit.aes = FALSE,
               aes(x=0,y=0, xend=PC1, yend=PC2, group=rowname),
               arrow = arrow(type = "closed",length=unit(4,'mm')),
               color= 'black') +
  geom_text(colour='black',size=4)+
  geom_text_repel(data= Veg.spe2.sc.ar, inherit.aes = FALSE,
                  aes(x= PC1,y=PC2,label=rowname))+
  labs(x = "PC1 30.5%", y = "PC2 21%")+ theme_minimal()

p
#setwd('graphs')
#ggsave(filename=paste('Veg_PCA_num.png'), width = 25, height = 15, dpi = 800, units = "cm", device='png')

plotly::ggplotly(p)

###Veg and Nute data####

Veg_Cover.num <- mutate_all(Veg_Cover.num, function(x) as.numeric(as.character(x)))
Veg_Nute_Join<-left_join(Nutrients.Sites,Veg_Cover.num, by="Site") #numeric site data
Veg_Nute_Join<-data.frame(Veg_Nute_Join, row.names=1)
str(Veg_Nute_Join)

Veg_Nute_Join_All<-left_join(Treat.w.nutes,Veg_Cover, by="Site")#all site data
colnames(Veg_Nute_Join_All)[3]<-'Interval'
colnames(Veg_Nute_Join_All)[2]<-'Severity'
Veg_Nute_Join_All_r<-data.frame(Veg_Nute_Join_All, row.names=1)

str(Veg_Nute_Join_All)

Veg_Nute.pca <- rda( Veg_Nute_Join~ Severity + Interval, data=Veg_Nute_Join_All_r, scale=TRUE)#RDAs explain 3.5%
Veg_Nute.pca <- rda( Veg_Nute_Join, data=Veg_Nute_Join_All_r, scale=TRUE)
###Making PCA#

plot(Veg_Nute.pca)    
summary(Veg_Nute.pca)

All.spe2.sc <- scores(Veg_Nute.pca, tidy=T)
All.spe2.sc.ar <- All.spe2.sc %>% 
  filter(score=='species')
All.spe2.sc.ar <- tibble::rownames_to_column(All.spe2.sc.ar)


All_Sites <- All.spe2.sc %>%
  filter(score=='sites') %>% 
  rename(Site=label) %>% 
  left_join(Site_Info)



p<-ggplot(All_Sites, aes(x=PC1, y=PC2, colour=Severity, label=Site )) + 
  geom_point(size=5.5)+ 
  geom_segment(data=All.spe2.sc.ar,inherit.aes = FALSE,
               aes(x=0,y=0, xend=PC1, yend=PC2, group=rowname),
               arrow = arrow(type = "closed",length=unit(4,'mm')),
               color= 'black') +
  geom_text(colour='black',size=4)+
  geom_text_repel(data= All.spe2.sc.ar, inherit.aes = FALSE, size=3,
                  aes(x= PC1,y=PC2,label=rowname))+
  labs(x = "PC1 25.1%", y = "PC2 22.7%")+ theme_minimal()

p
#setwd('graphs')
#ggsave(filename=paste('Nutes+Veg_PCA.png'), width = 25, height = 15, dpi = 800, units = "cm", device='png')


plotly::ggplotly(p)

