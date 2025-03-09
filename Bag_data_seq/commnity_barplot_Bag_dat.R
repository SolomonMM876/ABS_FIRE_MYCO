source('Bag_data_seq/Bag_ITS_data_prep.R')
library(ggplot2)



guild_perc<-dat_all_RA%>%
  group_by(guild2)%>%
  summarise(percent_guild2= (sum(count)/total_reads)*100)
guild_perc

taxa_w_guild<-all_tax%>%select(-guild2)%>%filter(!is.na(guild))

guild_comp<-all_tax%>%
  group_by(guild2)%>% 
  distinct()%>%
  summarise(OTUs= n())
guild_comp

Aphyla_comp<-all_tax%>%
  group_by(phylum)%>% 
  distinct()%>%
  summarise(OTUs= n())
phyla_comp

Phylum_comp<-dat_all_RA%>%
  group_by(phylum)%>%
  summarise(percent_phyla= (sum(count)/total_reads)*100)
Phylum_comp

min(dat_all_RA$reads_samp)
max(dat_all_RA$reads_samp)



library(patchwork)

#All guilds graphs
Bag_Interval_Comm<-dat_all_RA %>% 
  ggplot(aes(x=Fire.Interval, y=RA_total_interval, fill=guild2, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  labs(tag = "(a)") +
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18),
        plot.tag = element_text(size=25,hjust = 0.5),
        legend.position="none")+
         guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y='Proportion total sequencesfrom in-growth Bags', x= 'Fire Frequency') 
Bag_Interval_Comm


Bag_Severity_Comm<-dat_all_RA %>% 
  ggplot(aes(x=Fire.Severity, y=RA_total_severity, fill=guild2, text=species)) + # text aesthetic is for the ggplotly visualisation below
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
        legend.title = element_text(size = 18),
        legend.position="none")+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y= NULL, x= 'Fire Severity') 
Bag_Severity_Comm


Bag_Interval_Comm +Bag_Severity_Comm
#######################################################








#mycorrhizal explo_RA
Bag_Interval_explo<-dat_myco_RA %>% 
  ggplot(aes(x=Fire.Interval, y=RA_total_interval, fill=exploration_type, text=species)) + # text aesthetic is for the ggplotly visualisation below
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
  labs(y='Proportion EcM sequences from in-growth Bags', x= 'Fire Frequency') 
Bag_Interval_explo

Bag_Severity_explo<-dat_myco_RA %>% 
  ggplot(aes(x=Fire.Severity, y=RA_total_severity, fill=exploration_type, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18),
        legend.position="none")+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y= NULL, x= 'Fire Severity') 
Severity_explo

Bag_Interval_explo+Bag_Severity_explo
##################################################
#mycorrhizal genera

# Create a named vector for fill values
#this is extracted from community barplot Soil script
color_mapping <- setNames(genus_colors$fill_color, genus_colors$Genus)

Bag_Interval_Genus_myco<-dat_myco_RA %>% 
  ggplot(aes(x=Fire.Interval, y=RA_total_interval, fill=genus, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_fill_manual(values = color_mapping) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18),
        plot.tag = element_text(size=25,hjust = 0.5),
        legend.position="none")+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(tag = "(e)",y='Proportion EcM sequences from in-growth Bags', x= 'Fire Frequency') 
Bag_Interval_Genus_myco

Bag_Severity_Genus_myco<-dat_myco_RA %>% 
  ggplot(aes(x=Fire.Severity, y=RA_total_severity, fill=genus, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4)+
  scale_fill_manual(values = color_mapping) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y =element_blank(),
        legend.text = element_text(size = 12),  # Increase legend text size
        legend.title = element_text(size = 18),
        legend.position="none")+
    guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y= NULL, x= 'Fire Severity') 
Bag_Severity_Genus_myco

plotly::ggplotly(Bag_Severity_Genus_myco)

Bag_Interval_Genus_myco +Bag_Severity_Genus_myco
#################

#Assemble the Figure!
#(Bag_Interval_Comm +Bag_Severity_Comm )+(Bag_Interval_explo+Bag_Severity_explo)+ (Bag_Interval_Genus_myco +Bag_Severity_Genus_myco)


#### using both soil and bag

(Bag_Interval_Comm +Bag_Severity_Comm )|(Interval_Comm +Severity_Comm )->top
(Bag_Interval_explo+Bag_Severity_explo)|(Interval_explo+Severity_explo) ->middle
(Bag_Interval_Genus_myco +Bag_Severity_Genus_myco)|(Interval_Genus_myco +Severity_Genus_myco)->bottom
plot<-top/middle/bottom

ggsave(filename = "plots/Soil_Bag_community_plots.png", plot = plot, dpi=300, device = "png", width = 55, height = 60, units = "cm")