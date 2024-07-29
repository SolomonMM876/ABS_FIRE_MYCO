library(readxl)
library(dplyr)
library(ggplot2)
library(ggpmisc)

Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')


DW_DNA <- read_excel("Raw_data/DW_subsample_DNA_CNP.xlsx", 
                                     sheet = "Sheet2")%>%
  select(Tube_ID,DNA)%>%
  mutate(Tube_ID= as.character(Tube_ID))

library(readr)
qbit_DNA_Conc <- read_csv("Raw_data/Sequence/New_Sequencing/qbit_DNA_Conc_11_7_2024_ABS_Sites.csv")


DNA<-left_join(qbit_DNA_Conc%>%rename(Tube_ID= `Sample Name`),DW_DNA)%>%
  rename( Concentration_ng_uL= 'Sample Stock Concentration',
          mass_extracted= DNA)%>%
  left_join(Bag_Site%>%select(Site,Transect,Location,Tube_ID))%>%
  select(Site,Transect,Location,Tube_ID, mass_extracted,Concentration_ng_uL)

library(writexl)
write_xlsx(DNA,'Raw_data/Sequence/New_sequencing/DNA_ABS_STL_wt.xlsx' )

ggplot(DNA, aes(x = mass_extracted, y = Concentration_ng_uL)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE
  ) +
  labs(title = "DNA Weight vs Concentration",
       x = "Hyphal Biomass (mg)",
       y = "Concentration (ng/µL)") +
  theme_minimal()
