library(readxl)
library(dplyr)
library(ggplot2)
library(ggpmisc)

DW_DNA <- read_excel("Raw_data/DW_subsample_DNA_CNP.xlsx", 
                                     sheet = "Sheet2")%>%
  select(Tube_ID,DNA)%>%
  mutate(  Tube_ID= as.character(Tube_ID))

library(readr)
qbit_DNA_Conc <- read_csv("Raw_data/Sequence/qbit_DNA_Conc_11_7_2024_ABS_Sites.csv")


DNA<-left_join(qbit_DNA_Conc%>%rename(Tube_ID= `Sample Name`),DW_DNA)%>%
  rename( Concentration_ng_uL= 'Sample Stock Concentration')


ggplot(DNA, aes(x = DNA, y = Concentration_ng_uL)) +
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
