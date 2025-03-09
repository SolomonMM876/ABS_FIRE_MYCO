source('Soil_ITS_Scripts/Soil_ITS_data prep.R')
library(ggplot2)


trophic_mode<-dat_all_RA%>%
  group_by(trophicMode)%>%
  summarise(percent_troph= (sum(readcount)/total_reads)*100)
trophic_mode

taxa_w_guild<-tax%>%filter(!is.na(guild))

guild_comp<-tax%>%
  group_by(guild2)%>% 
  distinct()%>%
  summarise(OTUs= n())
guild_comp

phyla_comp<-tax%>%
  group_by(Phylum)%>% 
  distinct()%>%
  summarise(OTUs= n())
phyla_comp%>%select(Phylum)%>%pull()

Phylum_comp<-dat_all_RA%>%
  group_by(Phylum)%>%
  summarise(percent_phyla= (sum(readcount)/total_reads)*100)
Phylum_comp

library(patchwork)

#All guilds graphs
Interval_Comm<-dat_all_RA %>% 
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


Severity_Comm<-dat_all_RA %>% 
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


Interval_Comm +Severity_Comm
#######################################################








#mycorrhizal explo_RA
Interval_explo<-dat_myco_RA %>% 
  ggplot(aes(x=Interval, y=RA_total_interval, fill=exploration_type, text=Species)) + # text aesthetic is for the ggplotly visualisation below
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

Severity_explo<-dat_myco_RA %>% 
  ggplot(aes(x=Severity, y=RA_total_severity, fill=exploration_type, text=Species)) + # text aesthetic is for the ggplotly visualisation below
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

Interval_explo+Severity_explo
##################################################
#mycorrhizal genera

Interval_Genus_myco<-dat_myco_RA %>% 
  ggplot(aes(x=Interval, y=RA_total_interval, fill=Genus, text=Species)) + # text aesthetic is for the ggplotly visualisation below
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
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y='Proportion EcM sequences from Soil', x= 'Fire Frequency') 
Interval_Genus_myco

Severity_Genus_myco<-dat_myco_RA %>% 
  ggplot(aes(x=Severity, y=RA_total_severity, fill=Genus, text=Genus)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y =element_blank(),
        legend.text = element_text(size = 12),  # Increase legend text size
        legend.title = element_text(size = 18) )+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y= NULL, x= 'Fire Severity') 
Severity_Genus_myco


Interval_Genus_myco +Severity_Genus_myco

# Extract the ggplot build object
plot_build <- ggplot_build(Severity_Genus_myco)

# Extract the data
plot_data <- plot_build$data[[1]]

# Create a dataframe mapping Genus to Fill Colors
genus_colors <- unique(plot_data %>% select(fill, text))

# Rename columns for clarity
colnames(genus_colors) <- c("fill_color", "Genus")

#################

#Assemble the Figure!
(Interval_Comm +Severity_Comm )+(Interval_explo+Severity_explo)+ (Interval_Genus_myco +Severity_Genus_myco)

