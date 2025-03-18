library(readr)
library(tidyverse)


#read in community and bag data
#From Bag_ITS_data_prep
data<-read.csv('Processed_data/Bag_Seq_wide.csv')#bag data with selected meta and Myco communities


data<-data%>%filter(!Tube_ID %in% c(95,96))%>%#These samples are being removed because I pooled two locations from same Site/Tran during DNA extraction
  rename(readcount=myco_reads)%>%
  filter(!is.na(readcount)) # 6 samples are not included because they either did not contain mycorrhizal OTUs, did not have highly probable guild info or belonged to multiple guilds
  



#these data already have exploration types incorperated
#See BAG_ITS data prep
myco_dat<-read.csv('Processed_data/Bag_Seq_myco_dat.csv')


#read in genome size taxa df
tax_gs_bag<-read.csv('Processed_data/taxa_w_genome_size_bag_data.csv')
# 
 TrData<-myco_dat%>%
      left_join(tax_gs_bag%>%select(genus,mean_gs)%>%distinct())%>%
   distinct()




write.csv(data, file='HMSC_ABS/data/Bag_data.csv',row.names=FALSE)
write.csv(TrData, file='HMSC_ABS/data/Trait_Phylo_data.csv',row.names=FALSE)

###############
#CNP data prep
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')

Bag_data <- Bag_data %>%
  group_by(Site, Transect) %>%
  summarize(
    mean_log10_biomass_day = mean(log10_biomass_day, na.rm = TRUE),
    mean_Ortho_P_mg_kg = mean(Ortho_P_mg_kg, na.rm = TRUE),
    mean_Nitrate_mg_kg = mean(Nitrate_mg_kg, na.rm = TRUE),
    mean_Ammonia_mg_kg = mean(Ammonia_mg_kg, na.rm = TRUE),
    mean_pH = mean(pH, na.rm = TRUE)
  )

Stoich_Totals <- read.csv("Processed_data/CNP_clean.csv")%>%
  mutate(Site=as.numeric(Site),
         Transect=as.numeric(Transect))%>%
  rename(Carb_Hyph = Carbon, Nitrog_Hyph =Nitrogen, Phos_Hyph= Percent_Phos_)

Myco_abun<-read.csv('Processed_data/Myco_host_abundance.csv')

CNP_myco_comm<-read.csv('Processed_data/CNP_seq__myco_dat.csv')





CNP_data<-CNP_myco_comm%>%
group_by(Site,Transect, OTU) %>%
  #transform to wide format
  summarise(sequence_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(values_from = sequence_count, names_from = OTU,  values_fill = 0)%>%
  # Compute readcount as the sum of ITS columns
  mutate(readcount = rowSums(select(., starts_with("ITS")), na.rm = TRUE))%>%
  filter(!Site=='49')%>%
  left_join(Bag_data)%>%
  left_join(Stoich_Totals)%>%
  left_join(Myco_abun)%>%
  mutate( Site=as.factor(Site),
         Transect=as.factor(Transect))

TrData_CNP<-CNP_myco_comm %>%filter(!Site=='49')%>%
  select( kingdom:species,SH_species,OTU) %>% distinct()


TrData_CNP<-TrData_CNP%>%
  left_join(Fun_Traits, by = c('genus'='GENUS'))%>%
  select(kingdom:OTU,Ectomycorrhiza_exploration_type_template,Ectomycorrhiza_lineage_template)%>%
  rename(exploration_type=Ectomycorrhiza_exploration_type_template,
         Ecm_lineage=Ectomycorrhiza_lineage_template)%>%
  filter(!is.na(exploration_type))%>%
  left_join(tax_gs_bag%>%select(genus,mean_gs)
            %>%distinct())%>%filter(!is.na(mean_gs))
# I have to do this in order to make the phylo tree and run traits later, it is not a perfect solution, but it is the best I can think of


write.csv(CNP_data, file='HMSC_ABS/data/Bag_data_CNP.csv',row.names=FALSE)
write.csv(TrData_CNP, file='HMSC_ABS/data/Trait_Phylo_data_CNP.csv',row.names=FALSE)
