library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)
library(stringr)


DW_CNP <- read_excel("Raw_data/Stoich/DW_subsample_DNA_CNP.xlsx")
Bag_Site<-read.csv('Processed_data/All_Bag_Site_Info.csv')
#Nutrient resins
Resin_Nutrients<-read_excel('Processed_data/Resin_Nutrients.xlsx')

Manual_Calc_Total_P <- read_excel("Raw_data/Stoich/ABS_Total_P_2025.xlsx")

################## Second Round###########


###### Total_P ###
  
Total_P<-Manual_Calc_Total_P%>%
    rename( Sample =`Sample ID`,
           Percent_Phos_= `perc_P`,)%>%
    mutate( Sample = replace(Sample, row_number() %in% c(1, 5, 12), c(5.1, 10.2, 34.2)))%>%
   left_join(Bag_Site %>% dplyr::select(Site,Transect,Fire.Interval,Fire.Severity)%>% unique()%>%
                  mutate(Sample = paste(Site, Transect, sep = "."))
               , by = 'Sample')%>%
    filter(!str_detect(Sample, 'STAND|C C|MID|BLANK|blank|euc'))%>%
  dplyr::select(Site,Transect,Sample,Sample_mg, Fire.Interval, Fire.Severity,Percent_Phos_)%>%
  mutate(Site=as.factor(Site),
    Sample=as.factor(Sample),
    Site = case_when(
    !is.na(Site) ~ Site,  # Keep existing values in Site
    grepl("EUC", Sample) ~ "EUC",
    grepl("HFE", Sample) ~ "HFE",
    TRUE ~ Site           # Retain existing value if no condition matches
  ))%>%
  mutate(Fire.Interval= if_else(is.na(Fire.Interval),'Standards',Fire.Interval),
  Fire.Severity= if_else(is.na(Fire.Severity),'Standards',Fire.Severity),
  Stan_Group = ifelse(Fire.Interval == 'Standards', substr(Site, 1, 3), Fire.Interval),
  Percent_Phos_ = ifelse(Percent_Phos_ < 0 & Sample_mg > 1.1, 0.001/2, Percent_Phos_)) #lowest value on standard curve replaces neg value
  
  
  
p<-Total_P%>%
  #filter(Sample_mg>1)%>%
  ggplot(aes(x=Site, y=Percent_Phos_) )+ 
  geom_point(aes( color = Sample_mg < 1.1), size = 3,shape= 'square') +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  facet_grid(~Fire.Severity, scales = 'free_x')+
  ylab('Phosphorus %')+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))
p
plotly::ggplotly(p)
#CN
CNH <- read.csv("Processed_data/C_N_Stoich.csv")


CNH%>%
summarise(
  mean_Carbon = mean(Carbon, na.rm = TRUE),
  sd_Carbon = sd(Carbon, na.rm = TRUE),
  mean_Nitrogen = mean(Nitrogen, na.rm = TRUE),
  sd_Nitrogen = sd(Nitrogen, na.rm = TRUE),
  mean_C_N= mean(C_N),
  sd_C_N=sd(C_N)
)



interval_colors <- c("Long" = "darkred", "Short" = "orange")

CNH%>%
  filter(Sample_mg > .7)%>%
  filter(Fire.Interval!='Standards')%>%
  ggplot(aes(x=Fire.Interval) )+ 
  #geom_point(aes(y=C_N), size = 5, shape= 'triangle') +
  geom_boxplot(aes(fill=Fire.Interval, y= C_N))+
  scale_fill_manual(values = interval_colors) +     # Custom colors for Interval
  labs(x = "Fire Interval", y = "Hyphal Carbon:Nitrogen") +
  # annotate("text", x = 1.5, y = Inf, label = , hjust = 4, vjust = 1.1, size = 12)+
  theme_classic()+
  theme(axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        axis.line = element_line(size = 1.5),
        legend.position = 'none')
CNH%>%
  filter( !Fire.Interval=='Standards')%>%
ggplot(aes(x=as.factor(Site)) )+ 
  geom_point(aes(y=C_N), size = 5, shape= 'square') +
  geom_boxplot(aes(x=Fire.Severity, y= C_N))+
  facet_grid(~Fire.Severity, scales = 'free_x')+
  labs( y= ('C:N'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))



  min(CNH$C_N)


CNHP<-CNH%>%rename(Sample_mg_CN=Sample_mg)%>%mutate(Sample=as.factor(Sample))%>%
  full_join(Total_P%>%rename(Sample_mg_P=Sample_mg))

CNP_clean<-CNHP%>%
  # Sample 56.2 was excluded from analysis due to exceptionally low C and N values
  #Carbon= 40.39% Nitrogen= 1.50% Both values are  extremes of what was found in Zhang and Elser 2017 and lowest values found recorded for samples
  # Sampple is believed to contain 'other' non C/N material given the relatively low concentrations found
  filter(Sample!='56.2')%>%
  filter( !Fire.Interval=='Standards')%>%
  mutate(Percent_Phos_ = if_else(Percent_Phos_ < 0, NA_real_, Percent_Phos_))%>%
  select(-ID,-Stan_Group)%>%
  mutate(C_P = if_else(is.na(Percent_Phos_), NA_real_, Carbon / Percent_Phos_),
         N_P= Nitrogen/Percent_Phos_,
         C_N_P_ratio = paste(round(C_P, 0), round(N_P, 0), "1", sep = ":"))


CNP_clean%>%filter(!is.na(Percent_Phos_))

CNP_clean%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=C_P), size = 5, shape= 'square') +
 # geom_boxplot(aes(x=Fire.Severity, y= C_P))+
  labs( y= ('C:P'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

CNP_clean%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=C_N), size = 5, shape= 'square') +
  geom_boxplot(aes(x=Fire.Severity, y= C_N))+
  facet_grid(~Fire.Severity, scales = 'free_x')+
  labs( y= ('C:N'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))

CNP_clean%>%
  filter(!Site==12)%>%
  ggplot(aes(x=Site) )+ 
  geom_point(aes(y=N_P), size = 5, shape= 'square') +
  #geom_boxplot(aes(x=Fire.Interval, y= C_N_P))+
  #facet_grid(~Fire.Interval, scales = 'free_x')+
  labs( y= ('N:P'))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5))



write.csv(CNP_clean, "Processed_data/CNP_clean.csv",row.names = FALSE)

