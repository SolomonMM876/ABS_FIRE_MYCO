
library(ggplot2)
library(readxl)
library(dplyr)

Hyphae_CN <- read_excel("Raw_data/Stoich/Solomon Hyphae CN calc_11.07.24.xlsx", skip = 44)[-4,]
CN_Hyphae_ID <- read_excel("Raw_data/Stoich/CN_Hyphae.xlsx")
Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')


CNH<-Hyphae_CN%>%
  rename(Carbon=`Sample C calc %`,
         Nitrogen = `Sample N calc %`,
         Hydrogen= `Sample H calc %`,
         Sample_mg = `Sample mg`)%>%
  left_join(CN_Hyphae_ID%>%rename(ID=Well_ID))%>%mutate(Sample=as.character(Sample))%>%
  left_join(Bag_Site %>% select(Site,Transect, Fire.Interval, Fire.Severity)%>% unique()%>%
            mutate(Sample = paste(Site, Transect, sep = ".")))%>%
  select(Site,Transect,Sample,ID, Fire.Interval, Fire.Severity,Carbon, Nitrogen,Hydrogen,Sample_mg)%>%
  mutate(Site= if_else(is.na(Site),ID,Site),
         Fire.Interval= if_else(is.na(Fire.Interval),'Standards',Fire.Interval),
         Fire.Severity= if_else(is.na(Fire.Severity),'Standards',Fire.Severity),
         C_N = Carbon/Nitrogen)



CNH%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=Carbon, color = Sample_mg < .7), size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('% Carbon'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))


CNH%>%
  filter(!Sample==56.2)%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=Nitrogen, color = Sample_mg < .7), size = 3, shape= 'triangle') +
  geom_boxplot(aes(x=Fire.Interval, y= Nitrogen))+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('% Nitrogen'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

CNH%>%
  filter(!Sample==56.2)%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=C_N, color = Sample_mg < .7), size = 3, shape= 'triangle') +
  geom_boxplot(aes(x=Fire.Interval, y= C_N))+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('C:N'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))
