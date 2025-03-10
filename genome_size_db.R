
#LOAD LIBRARARIES

library(tidyverse)
library(vegan)
library(ggplot2)



#mycorrhizal taxa
myco_tax<-read.csv('Processed_data/Bag_dat_myco_tax.csv')

gs_dat<-read.csv('Processed_data/gs_dataset_mycorrhizae_Hiyang_3.3.25.csv')

gs_dat_genus<-gs_dat%>%filter(!is.na(genus))%>%group_by(genus)%>%summarise(mean_gs=mean(GS))

tax_gs_bag<-myco_tax%>%
left_join(gs_dat_genus)%>%filter(!is.na(mean_gs))


#Soil data
tax_myco_soil<-read.csv('Processed_data/Soil_ITS_myco_tax.csv')

tax_gs_soil<-tax_myco_soil%>%rename(genus=Genus)%>%
  left_join(gs_dat_genus)%>%filter(!is.na(mean_gs))




write.csv(tax_gs_bag, 'Processed_data/taxa_w_genome_size_bag_data.csv',row.names = FALSE)
write.csv(tax_gs_soil, 'Processed_data/taxa_w_genome_size_soil_data.csv',row.names = FALSE)

          