library(dplyr)
library(readr)
library(lubridate)
library(stringr)

NH4NO3_RESIN_LOT_1 <- read_csv("Raw_data/Nutes/24-07-10 SOLOMON NH4NO3 HCL RESIN LOT 1.csv")
NH4NO3_RESIN_LOT_2 <- read_csv("Raw_data/Nutes/24-07-12 SOLOMON NH4NO3 RESIN HCL LOT 2.csv")

NO3_NH4_Resin<-bind_rows(NH4NO3_RESIN_LOT_1,NH4NO3_RESIN_LOT_2)[,-c(10:13)] %>%
  rename(Test = `Test Name`, Sample_ID = `Sample ID`, Sample_Details = `Sample Details`) %>%
  mutate(Sample = case_when(
    !is.na(mdy(Sample_ID, quiet = TRUE)) ~ Sample_Details,
    TRUE ~ Sample_ID))%>%
  mutate(Sample = case_when(
    !is.na(dmy(Sample, quiet = TRUE)) ~ Sample_Details,
    TRUE ~ Sample))%>%
    mutate(Sample = case_when(
      is.na(Sample) ~ Sample_ID,
      TRUE ~ Sample))%>%
  select(-Time,-Sample_ID,-Sample_Details)


# Transform the data to wide format
wide_data <- NO3_NH4_Resin %>%
  pivot_wider(names_from = Test, values_from = c(Result,Absorbance,`Auto Dil`)) %>%
  rename(Ammonia= 'Result_Ammonia 2.0',
         Nitrate = 'Result_Nitrate 2',
         Ammonia_Dil = `Auto Dil_Ammonia 2.0`)%>%
  filter(!str_detect(Sample, 'STAND|C C|CC|NIRAJ'))%>%
  separate(Sample, into = c("Site", "Transect", "Location"), sep = "-", convert = TRUE)%>%
  mutate(Ammonia_mg_kg= Ammonia*7.5/1*Ammonia_Dil)#7.5mL used for 1 g of resin
  
  colnames(wide_data)

OPHOS_RESIN_HCL_LOT_1 <- read_csv("Raw_data/Nutes/24-07-13 SOLOMON OPHOS RESIN HCL LOT 1.csv")
OPHOS_RESIN_HCL_LOT_2 <- read_csv("Raw_data/Nutes/24-07-13 SOLOMON OPHOS RESIN HCL LOT 2.csv")

