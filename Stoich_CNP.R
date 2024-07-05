library(dplyr)
library(tidyr)
library(readr)
library(readxl)



DW_CNP <- read_excel("Raw_data/DW_subsample_DNA_CNP.xlsx")
Hyphal_P_1<-read_csv("Raw_data/Stoich/24-07-03 SOLOMON TKP LOT 1.csv")
Hyphal_P_2<-read_csv("Raw_data/Stoich/24-07-03 SOLOMON TKP LOT 2.csv")[-17,]

P_weights<-DW_CNP%>%
  select(Sample,P)%>%
  rename(Samp_mg=P)%>%
    mutate(Sample = replace(Sample, row_number() %in% c(1, 6, 8, 20), c(5.1, 8.2, 10.2, 34.2)))%>%
  left_join(
  #Hyphal_P_1
  Hyphal_P_1%>%
  rename(Sample='Sample ID',
         P_mg_L_1 = 'Adj Result',
         action = '...11',
         Final_Dil_1= 'Manual Dil')%>%
  mutate(Blank_1 =  mean(P_mg_L_1[Sample %in% c('2ND BLANK A','2ND BLANK B')],na.rm=T),
         Blank_Corrected_1= P_mg_L_1-Blank_1))%>%
  select(Sample, Blank_Corrected_1,action,Samp_mg,Final_Dil_1,P_mg_L_1)%>%
  #Hyphal_P_2
  left_join(Hyphal_P_2%>%
              select(`Sample ID`, `Adj Result`,`Manual Dil`, ...11 )%>%
  rename(Sample='Sample ID',
         P_mg_L_2 = 'Adj Result',,
         Final_Dil_2= 'Manual Dil')%>%
  mutate(Blank_2 =  mean(P_mg_L_2[Sample %in% c('BLANK')],na.rm=T),
         Blank_Corrected_2= P_mg_L_2-Blank_2))%>%
  mutate(Blank_Corrected = coalesce(Blank_Corrected_1, Blank_Corrected_2),
         action= coalesce(...11,action),
         Final_Dil= coalesce(Final_Dil_1,Final_Dil_2),
         P_mg_L = coalesce(P_mg_L_1,P_mg_L_2)) %>%
  select(Sample,Samp_mg,Blank_Corrected,Final_Dil,action,P_mg_L)%>%
  mutate( P_mg_Kg_Sol= (Blank_Corrected*Final_Dil*3/1000)/(Samp_mg/ 1000),
         Uncorrected_P_mg_Kg= (P_mg_L*Final_Dil*3/1000)/(Samp_mg/ 1000))

  library(writexl)

write_xlsx(P_weights, path = "Processed_data/concentration_P_hyph.xlsx")

