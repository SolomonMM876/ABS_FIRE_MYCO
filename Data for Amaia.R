library(dplyr)
library(tidyr)


dat <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_clean.tsv')

clean_dat <- dat %>%
  #remove transect pooled samples for CN
  filter(!is.na(Site))%>%
  #Samples 95 and 96 were pooled and have been separated below, but interpret with caution
  separate_rows(barcode, sep = "A") %>%
  mutate(barcode = ifelse(nchar(barcode) == 2, paste0("SMM", barcode), barcode))
  #SMM84 is both 11-1-16 and 11-2-16 

