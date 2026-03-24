
#LOAD LIBRARARIES

library(tidyverse)
library(vegan)
library(ggplot2)

#genome size db
gs_dat<-read.csv('Processed_data/gs_dataset_mycorrhizae_Hiyang_3.3.25.csv')

#mycorrhizal taxa
myco_tax<-read.csv('Processed_data/Bag_dat_myco_tax.csv')
myco_dat_bag<-read.csv( 'Processed_data/Bag_Seq_myco_dat.csv')




gs_dat_genus<-gs_dat%>%filter(!is.na(genus))%>%group_by(genus)%>%summarise(mean_gs=mean(GS))

gs_bag<-myco_dat_bag%>%
  left_join(gs_dat_genus)%>%filter(!is.na(mean_gs))%>%
  rename(readcount=count)%>%
  select(Site,Transect,readcount,genus,mean_gs)



#Soil data
tax_myco_soil<-read.csv('Processed_data/Soil_ITS_myco_tax.csv')
myco_dat_soil<-read.csv('Processed_data/Soil_ITS_myco_data.csv')


gs_soil<-myco_dat_soil%>%rename(genus=Genus)%>%
  left_join(gs_dat_genus)%>%filter(!is.na(mean_gs))%>%
  select(Site,Transect,readcount,genus,mean_gs)


hist(log10(gs_soil$mean_gs))
hist(log10(gs_bag$mean_gs))


# Combine datasets and label source
gs_combined <- bind_rows(
  gs_soil %>% mutate(source = "soil"),
  gs_bag %>% mutate(source = "bags")
)


ggplot(gs_combined, aes(x = source, y = mean_gs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # Boxplot without outliers
  geom_jitter(aes(size = readcount, color = source), width = 0.2, alpha = 0.7) + 
  scale_size_continuous(name = "Read Count") + 
  labs(title = "Comparison of Mean Genome Size", 
       x = "Source", 
       y = "Mean Genome Size") +
  theme_minimal()


# Compute weighted mean genome size per genus
# Compute relative abundance (RA) per genus within each source
gs_weighted <- gs_combined %>%
  group_by(source) %>%
  mutate(total_reads_source = sum(readcount, na.rm = TRUE)) %>%
  group_by(source, genus) %>%
  mutate(RA=readcount/total_reads_source,
         weight_gs= RA*mean_gs,
         log_weight_gs= log10(weight_gs))


t_test_result <- t.test(log_weight_gs ~ source, data = gs_weighted)

# Print test results
print(t_test_result)

library(ggpubr)


ggplot(gs_weighted, aes(x = source, y = 10^log_weight_gs, fill = source)) + 
  geom_violin( alpha = 0.5) +  # Boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.7, color = "black") +  # Jittered points
  scale_y_continuous(
    trans = "log10",  # Apply log10 scale to the y-axis for better visualization
    labels = scales::label_number(scale = 1)  # Format labels for easier interpretation
  ) +
  labs(title = paste('p=', round(t_test_result$p.value,2)),
       x = "Source",
       y = "Weighted Mean Genome Size") +
  theme_minimal()
