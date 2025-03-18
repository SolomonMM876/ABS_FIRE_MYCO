source('Soil_ITS_Scripts/Soil_ITS_data prep.R')
library(ggplot2)
library(patchwork)

######################################################
guild_perc<-dat_all_RA_soil%>%
  group_by(guild2)%>%
  summarise(percent_guild2= (sum(readcount)/first(total_reads))*100)
guild_perc


guild_comp<-tax%>%
  group_by(guild2)%>% 
  distinct()%>%
  summarise(OTUs= n())
guild_comp

phyla_comp<-tax%>%
  group_by(Phylum)%>% 
  distinct()%>%
  summarise(OTUs= n())
phyla_comp

Phylum_comp<-dat_all_RA_soil%>%
  group_by(Phylum)%>%
  summarise(percent_phyla= (sum(readcount)/first(total_reads))*100)
Phylum_comp

Phylum_comp<-dat_myco_RA_soil%>%
  group_by(Phylum,Interval)%>%
  summarise(percent_phyla= (sum(readcount)/first(interval_reads))*100)
Phylum_comp

Phylum_comp<-dat_myco_RA_soil%>%
  group_by(Phylum,Severity)%>%
  summarise(percent_phyla= (sum(readcount)/first(severity_reads))*100)
Phylum_comp

min(dat_all_RA_soil$reads_samp)
max(dat_all_RA_soil$reads_samp)

dat_all_RA_soil%>%select(Site,Transect,Severity,Interval)%>%group_by(Severity)%>%distinct()%>%summarise(n())

dat_all_RA_soil%>%select(Site,Transect)%>%distinct()%>%summarise(n())

#All guilds graphs
Interval_Comm<-dat_all_RA_soil %>% 
  ggplot(aes(x=Interval, y=RA_total_interval, fill=guild2, text=Species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  labs(tag = "(b)") +
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18),
        plot.tag = element_text(size=25,hjust = 0.5),
        legend.position="none")+
         guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y='Proportion total sequences from Soil', x= 'Fire Frequency') 
Interval_Comm


Severity_Comm<-dat_all_RA_soil %>% 
  ggplot(aes(x=Severity, y=RA_total_severity, fill=guild2, text=Species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  #scale_y_discrete(labels = NULL, breaks = NULL) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y= NULL, x= 'Fire Severity') 
Severity_Comm


#Interval_Comm +Severity_Comm
#######################################################








#mycorrhizal explo_RA
Interval_explo<-dat_explo_soil %>% 
  ggplot(aes(x=Interval, y=RA_explo_Interval, fill=exploration_type)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  labs(tag = "(d)") +
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18),
        plot.tag = element_text(size=25,hjust = 0.5),
        legend.position="none")+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y='Proportion EcM sequences from Soil', x= 'Fire Frequency') 
Interval_explo

Severity_explo<-dat_explo_soil%>% 
  ggplot(aes(x=Severity, y=RA_explo_Severity, fill=exploration_type)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y= NULL, x= 'Fire Severity') 
Severity_explo

#Interval_explo+Severity_explo
##################################################
#mycorrhizal genera

# Identify the top 10 most abundant genera
top_genera_soil <- dat_myco_RA_soil %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(RA_total_severity, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(Genus)

Interval_Genus_myco<-dat_myco_RA_soil %>% 
  mutate(genus_grouped = if_else(is.na(Genus) | !(Genus %in% top_genera_soil), "Other", Genus)) %>%
  ggplot(aes(x=Interval, y=RA_total_interval, fill=genus_grouped, text=Species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  labs(tag = "(f)") +
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18),
        plot.tag = element_text(size=25,hjust = 0.5),
        legend.position="none")+
  labs(y='Proportion EcM sequences from Soil', x= 'Fire Frequency') 
Interval_Genus_myco

Severity_Genus_myco<-dat_myco_RA_soil %>% 
  mutate(genus_grouped = if_else(is.na(Genus) | !(Genus %in% top_genera_soil), "Other", Genus)) %>%
  ggplot(aes(x=Severity, y=RA_total_severity, fill=genus_grouped, text=Genus)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y =element_blank(),
        legend.text = element_text(size = 12),  # Increase legend text size
        legend.title = element_text(size = 18) )+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 12))) +
  labs(y= NULL, x= 'Fire Severity') 
Severity_Genus_myco

#Interval_Genus_myco +Severity_Genus_myco

# Extract the ggplot build object
plot_build <- ggplot_build(Severity_Genus_myco)

# Extract the data
plot_data <- plot_build$data[[1]]

# Create a dataframe mapping Genus to Fill Colors
genus_colors <- unique(plot_data %>% select(fill, text))

# Rename columns for clarity
colnames(genus_colors) <- c("fill_color", "Genus")



