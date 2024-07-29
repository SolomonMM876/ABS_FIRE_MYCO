library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(stringr)


DW_CNP <- read_excel("Raw_data/DW_subsample_DNA_CNP.xlsx")
Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')
#Nutrient resins
Resin_Nutrients<-read_excel('Processed_data/Resin_Nutrients.xlsx')

Manual_Calc_Total_P <- read_excel("Raw_data/Stoich/Manual Calc_Total P.xlsx", 
                                  sheet = "2nd Run")

################## Second Round###########


###### Total_P ###
  
Total_P<-Manual_Calc_Total_P%>%
    rename( Sample =`Sample ID`,
           Percent_Phos_= `% P_ Sols calc`,)%>%
    mutate( Sample = replace(Sample, row_number() %in% c(1, 5, 12), c(5.1, 10.2, 34.2)))%>%
   left_join(Bag_Site %>% dplyr::select(Site,Transect,Fire.Interval,Fire.Severity)%>% unique()%>%
                  mutate(Sample = paste(Site, Transect, sep = "."))
               , by = 'Sample')%>%
    filter(!str_detect(Sample, 'STAND|C C|MID|BLANK'))%>%
  dplyr::select(Site,Transect,Sample,Sample_mg, Fire.Interval, Fire.Severity,Percent_Phos_)%>%
  mutate(Site= if_else(is.na(Site),Sample,Site),
  Fire.Interval= if_else(is.na(Fire.Interval),'Standards',Fire.Interval),
  Fire.Severity= if_else(is.na(Fire.Severity),'Standards',Fire.Severity),
  Stan_Group = ifelse(Fire.Interval == 'Standards', substr(Site, 1, 3), Fire.Interval))
  
  
p<-Total_P%>%
  #filter(Sample_mg>1)%>%
  ggplot(aes(x=Site, y=Percent_Phos_) )+ 
  geom_point(aes( color = Sample_mg < 1), size = 5,shape= 'square') +
  geom_boxplot(aes(x=Stan_Group))+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Interval, scales = 'free_x')+
  ylab('Phosphorus %')+
 # scale_y_continuous(breaks=seq(-1.5,2,by=.5))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

plotly::ggplotly(p)
#CN
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
summarise(
  mean_Carbon = mean(Carbon, na.rm = TRUE),
  sd_Carbon = sd(Carbon, na.rm = TRUE),
  mean_Nitrogen = mean(Nitrogen, na.rm = TRUE),
  sd_Nitrogen = sd(Nitrogen, na.rm = TRUE),
  mean_C_N= mean(C_N),
  sd_C_N=sd(C_N)
)


CNH%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=Carbon, color = Sample_mg < .7), size = 5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('% Carbon'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))


CNH%>%
  #filter(!Sample==56.2)%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=Nitrogen, color = Sample_mg < .7), size = 5, shape= 'triangle') +
  geom_boxplot(aes(x=Fire.Interval, y= Nitrogen))+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('% Nitrogen'))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

CNH%>%
  filter(Sample_mg > .7)%>%
  filter(Fire.Interval!='Standards')%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=C_N), size = 5, shape= 'triangle') +
  geom_boxplot(aes(x=Fire.Interval, y= C_N))+
 # scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('C:N'))+
  geom_hline(yintercept = 16, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25))+
  annotate("text", x = 6, y = Inf, label = "(F,P) =  (2.5317, 0.1461)", hjust = 1.1, vjust = 1.1, size = 5)
  



CNHP<-CNH%>%rename(Sample_mg_CN=Sample_mg)%>%
  full_join(Total_P%>%rename(Sample_mg_P=Sample_mg),
              by=c('Site','Transect','Fire.Interval','Fire.Severity','Sample'))%>%
  mutate(  C_P = Carbon/Percent_Phos_,
           C_N_P= Carbon/Nitrogen/Percent_Phos_)


CNHP%>%
  filter(Sample!='56.2')%>%
  filter( !Fire.Interval=='Standards')%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=C_P), size = 5, shape= 'square') +
  geom_boxplot(aes(x=Fire.Interval, y= C_P))+
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('C:P'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

CNHP%>%
  filter(Sample!='56.2')%>%
  filter( !Fire.Interval=='Standards')%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=C_N_P), size = 5, shape= 'square') +
  geom_boxplot(aes(x=Fire.Interval, y= C_N_P))+
  facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('C:N:P'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))



library(writexl)

write_xlsx(CNHP, path = "Raw_Data/Stoich/Stoich_Totals_Round_1.xlsx")








#Hyphal_P_1<-read_csv("Raw_data/Stoich/24-07-03 SOLOMON TKP LOT 1.csv")
#Hyphal_P_2<-read_csv("Raw_data/Stoich/24-07-03 SOLOMON TKP LOT 2.csv")[-17,]
#Hyphal_rerun<-read_csv("Raw_data/Stoich/24-07-08 SOLOMON TKP REPEAT.csv")
#First_Run<- read_excel("Processed_data/concentration_P_hyph.xlsx")

#Im not using this data, but dont want to erase the script



# #First Run
# P_weights<-DW_CNP%>%
#   select(Sample,P)%>%
#   rename(Samp_mg=P)%>%
#     mutate(Sample = replace(Sample, row_number() %in% c(1, 6, 8, 20), c(5.1, 8.2, 10.2, 34.2)))%>%
#   left_join(
#   #Hyphal_P_1
#   Hyphal_P_1%>%
#   rename(Sample='Sample ID',
#          P_mg_L_1 = 'Adj Result',
#          action = '...11',
#          Final_Dil_1= 'Manual Dil')%>%
#   mutate(Blank_1 =  mean(P_mg_L_1[Sample %in% c('2ND BLANK A','2ND BLANK B')],na.rm=T),
#          Blank_Corrected_1= P_mg_L_1-Blank_1))%>%
#   select(Sample, Blank_Corrected_1,action,Samp_mg,Final_Dil_1,P_mg_L_1)%>%
#   #Hyphal_P_2
#   left_join(Hyphal_P_2%>%
#               select(`Sample ID`, `Adj Result`,`Manual Dil`, ...11 )%>%
#   rename(Sample='Sample ID',
#          P_mg_L_2 = 'Adj Result',,
#          Final_Dil_2= 'Manual Dil')%>%
#   mutate(Blank_2 =  mean(P_mg_L_2[Sample %in% c('BLANK')],na.rm=T),
#          Blank_Corrected_2= P_mg_L_2-Blank_2))%>%
#   mutate(Blank_Corrected = coalesce(Blank_Corrected_1, Blank_Corrected_2),
#          action= coalesce(...11,action),
#          Final_Dil= coalesce(Final_Dil_1,Final_Dil_2),
#          P_mg_L = coalesce(P_mg_L_1,P_mg_L_2)) %>%
#   select(Sample,Samp_mg,Blank_Corrected,Final_Dil,action,P_mg_L)%>%
#   mutate( P_mg_Kg_Sol= (Blank_Corrected*Final_Dil*3/1000)/(Samp_mg/ 1000),
#          Uncorrected_P_mg_Kg= (P_mg_L*Final_Dil*3/1000)/(Samp_mg/ 1000))
# 


# mutate( P_mg_Kg_Sol= (Blank_Corrected*Final_Dil*3/1000)/(Samp_mg/ 1000),
#         Uncorrected_P_mg_Kg= (P_mg_L*Final_Dil*3/1000)/(Samp_mg/ 1000))

