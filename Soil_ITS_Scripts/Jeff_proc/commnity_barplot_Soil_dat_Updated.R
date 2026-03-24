source('Soil_ITS_Scripts/Jeff_Proc/Soil_ITS_data_prep_Jeff_up.R')
library(ggplot2)


dat_all_RA_soil%>%  filter(str_detect(guild, "mycorrhizal") | str_detect(guild, "Mycorrhizal")) %>%
  select(OTU)%>%distinct()


######################################################
guild_perc<-dat_all_RA_soil%>%
  group_by(guild2)%>%
  summarise(percent_guild2= (sum(count)/first(total_reads))*100)
guild_perc


guild_comp<-all_tax_soil%>%
  group_by(guild2)%>% 
  distinct()%>%
  summarise(OTUs= n())
guild_comp

phyla_comp<-all_tax_soil%>%
  group_by(phylum)%>% 
  distinct()%>%
  summarise(OTUs= n())
phyla_comp

Phylum_comp<-dat_myco_RA_soil%>%
  group_by(phylum)%>%
  summarise(percent_phyla= (sum(count)/first(myco_reads))*100)
Phylum_comp

Phylum_comp<-dat_myco_RA_soil%>%
  group_by(phylum,Fire.Interval)%>%
  summarise(percent_phyla= (sum(count)/first(interval_reads))*100)
Phylum_comp

Phylum_comp<-dat_myco_RA_soil%>%
  group_by(phylum,Fire.Severity)%>%
  summarise(percent_phyla= (sum(count)/first(severity_reads))*100)
Phylum_comp

min(dat_all_RA_soil$reads_samp)
max(dat_all_RA_soil$reads_samp)

dat_all_RA_soil%>%select(Site,Transect,Fire.Severity,Fire.Interval)%>%group_by(Fire.Severity)%>%distinct()%>%summarise(n())

dat_all_RA_soil%>%select(Site,Transect)%>%distinct()%>%summarise(n())

#All guilds graphs
Interval_Comm<-dat_all_RA_soil %>% 
  ggplot(aes(x=Fire.Interval, y=RA_total_interval, fill=guild2, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  scale_fill_viridis_d(option = "D", name = "Guild") +  # legend title here
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
  labs(y='Proportion total sequences', x= 'Fire Frequency') 
Interval_Comm


Severity_Comm<-dat_all_RA_soil %>% 
  ggplot(aes(x=Fire.Severity, y=RA_total_severity, fill=guild2, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  scale_fill_viridis_d(option = "D", name = "Guild") +  # legend title here
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18))+
  guides(fill = guide_legend(override.aes = list(shape = 14, size = 15))) +
  labs(y= NULL, x= 'Fire Severity') 
Severity_Comm


#Interval_Comm +Severity_Comm
#######################################################








#mycorrhizal explo_RA
# Interval_explo<-dat_explo_soil %>% 
#   ggplot(aes(x=Fire.Interval, y=RA_explo_Interval, fill=exploration_type)) + # text aesthetic is for the ggplotly visualisation below
#   geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
#   scale_x_discrete(drop=FALSE) + 
#   theme_classic()+
#   labs(tag = "(d)") +
#   theme(axis.text.x = element_text(hjust = 0.5,size=20),
#         axis.text.y = element_text(size=20),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_text(size=25),
#         legend.text = element_text(size = 15),  # Increase legend text size
#         legend.title = element_text(size = 18),
#         plot.tag = element_text(size=25,hjust = 0.5),
#         legend.position="none")+
#   guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
#   labs(y='Proportion EcM sequences from Soil', x= 'Fire Frequency') 
# Interval_explo
# 
# Severity_explo<-dat_explo_soil%>% 
#   ggplot(aes(x=Fire.Severity, y=RA_explo_Severity, fill=exploration_type)) + # text aesthetic is for the ggplotly visualisation below
#   geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
#   scale_x_discrete(drop=FALSE) + 
#   theme_classic()+
#   theme(axis.text.x = element_text(hjust = 0.5,size=20),
#         axis.text.y = element_blank(),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_blank(),
#         legend.text = element_text(size = 15),  # Increase legend text size
#         legend.title = element_text(size = 18))+
#   guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
#   labs(y= NULL, x= 'Fire Severity') 
# Severity_explo

#Fire.Interval_explo+Severity_explo
##################################################
#mycorrhizal genera

# Identify the top 10 most abundant genera
top_genera_soil <- dat_myco_RA_soil %>%
  group_by(genus) %>%
  summarise(total_abundance = sum(RA_total_severity, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(genus)
# Compute ordered levels for genus_grouped
genus_order <- dat_myco_RA_soil %>%
  mutate(genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus)) %>%
  group_by(genus_grouped) %>%
  summarise(total_abundance = sum(RA_total_severity, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%  # ascending so the most abundant ends up on top
  pull(genus_grouped)

genus_order <- c(setdiff(genus_order, "Other"), "Other")

# Make genus_grouped a factor with the correct order
dat_myco_RA_soil <- dat_myco_RA_soil %>%
  mutate(genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus),
         genus_grouped = factor(genus_grouped, levels = genus_order))

Interval_Genus_myco<-dat_myco_RA_soil %>% 
  ggplot(aes(x=Fire.Interval, y=RA_total_interval, fill=genus_grouped, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_fill_viridis_d(option = "D", name = "Genus") +  # legend title here
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
  labs(y='Proportion Mycorrhizal sequences', x= 'Fire Frequency') 
Interval_Genus_myco

Severity_Genus_myco<-dat_myco_RA_soil %>% 
  ggplot(aes(x=Fire.Severity, y=RA_total_severity, fill=genus_grouped, text=genus)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  scale_fill_viridis_d(option = "D", name = "Genus") +  # legend title here
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=25),
        axis.title.y =element_blank(),
        legend.text = element_text(size = 17),  # Increase legend text size
        legend.title = element_text(size = 18) )+
  guides(fill = guide_legend(override.aes = list(shape = 14, size = 15))) +
  labs(y= NULL, x= 'Fire Severity') 
Severity_Genus_myco

#Interval_Genus_myco +Severity_Genus_myco



