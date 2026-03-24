library(readxl)
Hyphal_Pixels <- read_excel("Raw_data/ABS_Hyphal_Data.xlsx")

Hyphal_um <- Hyphal_Pixels %>%
  #100um/884.366 pixels = 1um/8.84366
  mutate(    across(Line1:`Line 20`, ~ as.numeric(.) * (1/8.84366)))

    10*1/8.84366
    
    