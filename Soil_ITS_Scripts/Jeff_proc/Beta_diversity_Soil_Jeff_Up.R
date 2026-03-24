#LOAD LIBRARARIES
library(forcats)
library(tidyr)
library(dplyr)
library(vegan)
library(indicspecies)
library(tibble) 
library(readxl)
library(ggplot2)



#From ITS prep script
Soil_Seq_wide<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_Seq_wide.csv')
myco_tax<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_dat_myco_tax.csv')


# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from one or more groups of veg class
multipatt(Soil_Seq_wide %>% select(starts_with('ITSall')), # first argument is the community table, select only those columns
          Soil_Seq_wide$Fire.Severity) -> res.severity

multipatt(Soil_Seq_wide %>% select(starts_with('ITSall')), # first argument is the community table, select only those columns
          Soil_Seq_wide$Fire.Interval) -> res.interval

summary(res.interval)
summary(res.severity)

total_myco_reads<-Soil_Seq_wide%>% select(starts_with('ITSall'))%>% 
  pivot_longer(cols=starts_with('ITSall'), names_to='OTU', 
               values_to='count')%>%left_join(myco_tax)%>%
  summarise(sum(count))%>%
  pull()

total_myco_reads

# to visualise differences in taxonomic composition, using output from multipatt()
# first prepare the data in the res object so that it can be joined with taxonomic info
# output is only those otus significant for one or more pairs
out.interval <- res.interval[['sign']] %>% #what does this do?
  filter(p.value <= 0.05) %>% 
  rownames_to_column('OTU') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(site_code = gsub('^s.', '', group))

out.severity <- res.severity[['sign']] %>% #what does this do?
  filter(p.value <= 0.05) %>% 
  rownames_to_column('OTU') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(site_code = gsub('^s.', '', group))
#then the relevant community data, 
# and reorder the otu levels by decreasing abundance
out.severity <- left_join(out.severity, Soil_Seq_wide %>% 
              select(Site,Transect, Fire.Severity, starts_with('ITSall'))%>% 
              pivot_longer(cols = starts_with('ITSall'), names_to = 'OTU', values_to = 'count'))%>%
                group_by(OTU,Site,Transect,Fire.Severity,count)%>%
                summarize(
                  count=sum(count))%>%
  mutate(OTU = fct_reorder(OTU, count, max), 
         Fire.Severity = as_factor(Fire.Severity))%>%
  left_join(myco_tax)


out.interval <- left_join(out.interval, Soil_Seq_wide %>% 
              select(Site,Transect,Fire.Interval,starts_with('ITSall'))%>% 
  pivot_longer(cols = starts_with('ITSall'), names_to = 'OTU', values_to = 'count'))%>%
  mutate(OTU = fct_reorder(OTU, count, max), 
         Interval = as_factor(Fire.Interval))%>%
  left_join(myco_tax)%>%
  filter(count>0)

# then join with the taxonomy table
tax.interval <- left_join(out.interval, myco_tax)
tax.interval

out.interval.filter<-out.interval%>%
  filter(count>0)%>%
  mutate( rel_abun= (count/total_myco_reads)*100)


# finally produce the barplot
out.severity.filter<-out.severity%>%
  filter(count>0)%>%
  mutate( rel_abun= (count/total_myco_reads)*100)

p1<-out.severity.filter %>% 
  ggplot(aes(x=Fire.Severity, y=rel_abun, fill=genus, text=OTU)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity', position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  #scale_fill_manual(values = custom_palette) +  #palette.pals()
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18) )+ # Increase legend title size)+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 15 ))) +
  labs(y='Relative Abundance %', x= 'Fire Severity') 
p1

p2<-out.interval.filter %>% 
  ggplot(aes(x=Fire.Interval, y=rel_abun, fill=genus, text=OTU)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity', position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  #scale_fill_manual(values = custom_palette) +  #palette.pals()
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18) )+ # Increase legend title size)+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 15 ))) +
  labs(y='Relative Abundance %', x= 'Fire Frquency') 
p2

  
library(patchwork)
p1+p2

# use this next lines to interactively get OTU IDs
plotly::ggplotly(p1)
#######################


###next analysis#######

# next analysis - permanova

# first remove samples that have no mycorrhizal reads because they cause errors below
Soil_Seq_wide<-Soil_Seq_wide%>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

# extract the community table, save as a new object
mat_myco<- Soil_Seq_wide %>% select(starts_with("ITSall"))



temp<-adonis2(mat_myco~ Fire.Severity*Fire.Interval , data=Soil_Seq_wide, distance='robust.aitchison', by= 'margin')

temp<-temp%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>% 
  mutate(Sample_Type='Soil')
temp
write.csv(temp,'Tables/Soil_permanova.csv', row.names = FALSE)
table(Soil_Seq_wide$Fire.Interval)
table(Soil_Seq_wide$Fire.Severity)

# adonis2(mat_ecm_gs~ Fire.Severity+Fire.Interval + Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
#           log10_biomass_day+myco_host_total+gs_mean
#         ,data=Soil_Seq_wide, distance='robust.aitchison', add=TRUE)


cap.all <- capscale(mat_myco~ Fire.Interval * Fire.Severity 
                    , data=Soil_Seq_wide, distance='robust.aitchison', add=TRUE)
anova(cap.all)
summary(cap.all)
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
Cap1_aov
cap.all
plot(cap.all)
proportions<-round(cap.all$CCA$eig/cap.all$tot.chi *100, 1) # proportion of variation associated with each axis

#cap_test <- capscale(mat_myco ~ 1 , data=Soil_Seq_wide, distance='robust.aitchison', add=TRUE)

#ordistep(cap_test, formula(cap.all), direction='forward')



# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap.all, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site<- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')
scrs_biplot <- scrs %>% filter(score=='biplot')


#Trying different scaling
# scrs_scaled <- scores(cap.all, scaling = "sites", tidy=TRUE)
# scrs_scaled<- scores(cap.all, scaling = "species", tidy=TRUE)
# scrs_scaled <- scores(cap.all, scaling = "symmetric", tidy=TRUE)
# 
# 
# scrs_biplot <- scrs_scaled %>% filter(score=='biplot')
# scrs_site <- scrs_scaled %>% filter(score=='sites')






interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_colors <- c("Hgih" = "blue", "Low" = "green")

library(ggrepel)

# first plot - site scores along with centroids for each group
cbind(Soil_Seq_wide, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2)) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_point(aes( colour= Fire.Interval,shape= Fire.Severity), size=10, stroke = 3)+ 
  stat_ellipse(aes(color = Fire.Interval), level = 0.95, size = 2, linetype = 1) + # Add confidence ellipses
  #geom_text(aes( label = Site), color= 'black', size=4)+
  #geom_text(data = scrs_cent, aes(label = label), size = 2) + 
  scale_shape_manual(values = c(19,1))+
  scale_colour_manual(values = interval_colors) +     # Custom colors for Interval
  labs( x=  paste0("CAP1 (", proportions[1], "%)"), y=  paste0("CAP2 (", proportions[2], "%)"),
              colour = "Fire Frequency",  # Rename legend for color (Fire.Interval)
              shape = "Fire Severity"     # Rename legend for shape (Fire.Severity)
        )+
  # geom_segment(data=scrs_biplot%>%
  #                filter(abs(CAP1) > 0.3 | abs(CAP2) > 0.3),
  #              inherit.aes = FALSE,
  #              aes(x=0,y=0, xend=CAP1, yend=CAP2, group=label),
  #              arrow = arrow(type = "closed",length=unit(3,'mm')),
  #              color= 'black') +
  # geom_text_repel(data=scrs_biplot%>%
  #                   filter(abs(CAP1) > 0.3 | abs(CAP2) > 0.3),
  #                 aes(x=CAP1, y=CAP2, label=label),
  #                 colour='black',size=6)+
  xlim(c(min(scrs_site[, 'CAP1']-1), max(scrs_site[, 'CAP1'])+.2)) + 
  ylim(c(min(scrs_site[, 'CAP2']), max(scrs_site[, 'CAP2'])+1)) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 36, color = 'black'),
    axis.text.y = element_text(size = 36, color = 'black'),
    axis.title.y = element_text(size = 36, color = 'black'),
    axis.title.x = element_text(size = 36),
    legend.position = "top",
    legend.text = element_text(size = 30),        # Larger legend text
    legend.title = element_text(size = 30)        # Optional: larger legend title
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 8)),  # Adjust shape size
    colour = guide_legend(override.aes = list(size = 8))  # Adjust color legend size
  )  ->p2

p2

ggsave(filename = "plots/30_Soil_PCoA.png", plot = p2, dpi=300, device = "png", width = 65, height = 45, units = "cm")

