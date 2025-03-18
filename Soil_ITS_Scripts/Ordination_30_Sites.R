source('Soil_ITS_Scripts/Soil_ITS_data prep.R')

#LOAD LIBRARARIES
library(tidyverse)
library(vegan)
library(indicspecies)

#select only mycorrhizal OTUs from sites we previously selected
mat_myco <-  filtered_mat %>% as.data.frame() %>% 
  select(all_of(myco_otus))%>%# just ecto OTUs
  filter(rownames(.) %in% Site_Sample_IDs)  # Keep only the desired samples


dat_ecm_30_site <- left_join( mat_myco%>%
                        rownames_to_column('sample_ID'), Blast_ID) 


# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from fire interval
res_Interval<-multipatt(dat_ecm_30_site%>% select(ends_with('.09FU')), # first argument is the community table, select only those columns
                        dat_ecm_30_site$Interval) 
summary(res_Interval)


res_Severity<-multipatt(dat_ecm_30_site%>% select(ends_with('.09FU')), # first argument is the community table, select only those columns
                        dat_ecm_30_site$Severity) 
summary(res_Severity)


# to visualise differences in taxonomic composition, using output from multipatt()
# first prepare the data in the res object so that it can be joined with taxonomic info
# Interval
out_Interval <- res_Interval[['sign']] %>% 
  filter(p.value <= 0.05) %>% 
  rownames_to_column('SH_ID') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(Interval = gsub('^s.', '', group))

#Severity
out_Severity <- res_Severity[['sign']] %>% #what does this do?
  filter(p.value <= 0.05) %>% 
  rownames_to_column('SH_ID') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(Severity = gsub('^s.', '', group))
# then join with the taxonomy table, then the relevant community data, 
# and reorder the otu levels by decreasing abundance
out_Interval <- left_join(out_Interval, tax) %>% 
  left_join(dat_ecm_30_site %>% 
              select(Site,Transect, sample,sample_ID, Interval, ends_with('.09FU')) %>% 
              pivot_longer(cols=ends_with('.09FU'), names_to='SH_ID', 
                           values_to='count')) %>% 
  mutate(SH_ID = fct_reorder(SH_ID, count, max), 
         Interval = as_factor(Interval))

out_Severity <- left_join(out_Severity, tax) %>% 
  left_join(dat_ecm_30_site %>% 
              select(Site,Transect, sample,sample_ID, Severity, ends_with('.09FU')) %>% 
              pivot_longer(cols=ends_with('.09FU'), names_to='SH_ID', 
                           values_to='count')) %>% 
  mutate(SH_ID = fct_reorder(SH_ID, count, max), 
         Severity = as_factor(Severity))

library(RColorBrewer)


# Generate a custom color palette by combining multiple RColorBrewer palettes
custom_palette <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(9, "Set1"))

library(Polychrome)
set.seed(723451) # for reproducibility
custom_palette <- createPalette(17, c("#ff0000"), M=1000)
swatch(custom_palette)

#calculate relative abundance
p<-dat_ecm_30_site%>% select(ends_with('.09FU'))%>% 
  pivot_longer(cols=ends_with('.09FU'), names_to='SH_ID', 
               values_to='count')%>%left_join(tax)%>%
  summarise(sum(count))%>%
  pull()

p


# finally produce the barplot
Interval_Indicator<-out_Interval %>% 
  ggplot(aes(x=Interval, y=(count/p)*100, fill=Genus, text=SH_ID)) + # text aesthetic is for the ggplotly visualization below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  scale_x_discrete(drop=FALSE) + 
  theme_classic()+
 # scale_fill_manual(values = custom_palette) + 
  theme(axis.text.x = element_text(hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25),
        legend.text = element_text(size = 15),  # Increase legend text size
        legend.title = element_text(size = 18) )+
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 8))) +
  labs(y='Relative Abundance %', x= 'Fire Interval') 

Interval_Indicator

Severity_Indicator<-out_Severity %>% 
  ggplot(aes(x=Severity, y=(count/p)*100, fill=Genus, text=SH_ID)) + # text aesthetic is for the ggplotly visualisation below
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

Severity_Indicator

library(patchwork)
Interval_Indicator+Severity_Indicator
# use this next lines to interactively get OTU IDs
plotly::ggplotly(Interval_Indicator)
#######################




####next analysis#######

# next analysis - permanova
# extract the community table, save as a new object

mat_ecm <- dat_ecm_30_site %>% select(ends_with('.09FU'))



adonis2(mat_ecm ~ Severity+Interval , data=dat_ecm_30_site, distance='robust.aitchison', add=TRUE)

table(dat_ecm_30_site$Interval)
table(dat_ecm_30_site$Severity)



cap.all <- capscale(mat_ecm~ Interval+Severity  , data=dat_ecm_30_site, distance='robust.aitchison', add=TRUE)



# a constrained analysis of principal coordinates using a different distance index - result is quite good
# cap1 <- capscale(mat_ecm ~ Interval+Severity +Veg_Class +
#                    Tree.Basal.Area_m2 + Herb.Cover_0.50cm_perc + perc_myco_host_freq+
#                    NH4 + NO3 + Total.P
#                  , data=dat_ecm_30_site, distance='robust.aitchison', add=TRUE)
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
plot(cap.all)
cap.all # summary of inertia
proportions<-round(cap.all$CCA$eig/cap.all$tot.chi *100, 1) # proportion of variation associated with each axis
proportions
anova(cap.all) # statistical significance of the constraint

#cap_test <- capscale(mat_ecm ~ 1 , data=dat_ecm_30_site, distance='robust.aitchison', add=TRUE)
#ordistep(cap_test, formula(cap1), direction='forward')



# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap.all, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')%>%
  filter(abs(CAP1) > 0.6 | abs(CAP2) > 0.6)%>%
  rename(SH_ID=label)%>%
  left_join(tax)
scrs_site <- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')
scrs_biplot <- scrs %>% filter(score=='biplot')


interval_colors <- c("Long" = "darkred", "Short" = "orange")
#scrs_cent$label <- c("Interval (Long)", "Interval (Short)", "Severity (High)", "Severity (Low)")


# first plot - site scores along with centroids for each group
p2 <- cbind(dat_ecm_30_site, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2)) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_point(aes(colour= Interval, shape= Severity), size=8, stroke = 3) + 
  stat_ellipse(aes(color = Interval), level = 0.95, size = 2, linetype = 1) + # Add confidence ellipses
  scale_shape_manual(values = c(19,1))+
  scale_colour_manual(values = interval_colors) +     
  labs(x= paste0("CAP1 (", proportions[1], "%)"), y= paste0("CAP2 (", proportions[2], "%)")) +
  geom_segment(data=scrs_biplot %>% filter(abs(CAP1) > 0.1 | abs(CAP2) > 0.1),
               inherit.aes = FALSE,
               aes(x=0, y=0, xend=CAP1, yend=CAP2, group=label),
               arrow = arrow(type = "closed", length=unit(3, 'mm')),
               color= 'black') +
  geom_text(data=scrs_biplot %>% filter(abs(CAP1) > 0.1 | abs(CAP2) > 0.1),
                  inherit.aes = FALSE,
                  aes(x=CAP1, y=CAP2, label=label),
                  colour='black', size=10, fontface="bold") +
  xlim(c(min(scrs_site[, 'CAP1']), max(scrs_site[, 'CAP1'])+.5)) + 
  ylim(c(min(scrs_site[, 'CAP2'])-0.5, max(scrs_site[, 'CAP2'])+2)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5, size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25)) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 20)),
         shape = guide_legend(override.aes = list(color = "black", size = 20))) +
  theme(legend.position='top')
  
p2


#left_join(Bag_Site %>% select(Site,Transect, Site_Pair)%>% unique(), by = c("Site","Transect"))
  

# plot side-by-side using the patchwork package
#library(patchwork)
#p1 + p2


# top ten otus associated with the top of the ordination, presumably '???' samples
otus<-scrs_spp%>%
  arrange(desc(CAP2)) %>% 
  head(10)
# taxonomic information for those otus
tax %>% 
  filter(SH_ID %in% otus$label)


# still tidying - plotting turnover in space
pco1 <- capscale(mat_ecm ~ 1, data=dat_ecm_30_site, distance='robust.aitchison', add=TRUE)
scrs_site <- scores(pco1, display='sites') # TIDY
cbind(dat_ecm_30_site, scrs_site) %>% 
  ggplot(aes(x=Longitude, y=Latitude, colour=MDS1)) + 
  geom_point()


# including spatial patterns in analyses of community composition
# one way to do this is using principle coordinates of neighbour matrices
xy <- dist(dat_ecm_30_site[, c('Longitude', 'Latitude')]) # create distance matrix of longs and lats
xy.pcnm <- pcnm(xy)$vectors %>% as.data.frame() # export the result into a dataframe

# combine the PCNM results with the rest of the data  
temp <- cbind(dat_ecm_30_site, xy.pcnm)

# what do the PCNMs represent? A plotting example
ggplot(temp, aes(x=Longitude, y=Latitude, colour=PCNM1)) + 
  geom_point()

# what variation is explained by these spatial variables
cap.sp <- capscale(mat_ecm ~ ., data=temp %>% 
                     select('Longitude', 'Latitude', starts_with('PCNM')))
cap.sp
# proportion of variation associated with each axis
cap.sp$CCA$eig/cap.sp$tot.chi
# explanatory value
anova(cap.sp)

# which individual spatial variables to include?
cap.0 <- capscale(mat_ecm ~ 1, data=temp) # intercept-only, starting analysis
# uncomment this next line to run - takes a long time
#ordistep(cap.0, formula(cap.sp), direction='forward')


# should we include all climate variables in our analysis
# check for variance inflation - values higher than ~ 10 are unlikely to explain unique variation
#cap.cl <- capscale(mat_ecm ~ Annual_Prec + elev, data=dat_ecm_30_site)
#vif.cca(cap.cl)

# variation partitioning - here to three groups of variables
vp <- varpart(vegdist(mat_ecm, distance='robust.aitchison'), 
              ~Interval+Severity, # Most sig factors tested, though Interval is actually sig = X1
              ~   NH4 + NO3 + Bray.P+Total.P , # climate = X2
              ~ Tree.Basal.Area_m2 + Herb.Cover_0.50cm_perc + perc_myco_host_freq, # veg = X3
              ~ Site # Effect of site = X4
            , data=temp)
vp
plot(vp, Xnames=c('Regime', 'Nutrients','Veg','Site'))

# test unique variation explained by partitions using partial CAPs
# 1 - associated with Fire REgime
anova(capscale(mat_ecm ~ Interval + Total.P + 
                 Condition(Annual_Prec + elev +
                             PCNM7 + PCNM1 + PCNM8 + PCNM4 ), data=temp, 
               distance='robust.aitchison'))
# 2 - associated with climate
anova(capscale(mat_ecm ~ Annual_Prec+ elev +
                 Condition(Interval + Total.P + 
                             PCNM7 + PCNM1 + PCNM8 + PCNM4), data=temp, 
               distance='robust.aitchison'))
# 3 - associated with space
anova(capscale(mat_ecm ~ PCNM7 + PCNM1 + PCNM8 + PCNM4 +
                 Condition(Annual_Prec+ elev + Interval + Total.P), data=temp, 
               distance='robust.aitchison'))




