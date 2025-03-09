#LOAD LIBRARARIES
library(forcats)
library(tidyr)
library(dplyr)
library(vegan)
library(indicspecies)
library(tibble) 
library(readxl)
library(ggplot2)



Bag_Seq_wide<-read.csv('Processed_data/Bag_Seq_wide.csv')
Bag_data<-read.csv('Processed_data/Updated_Bag_data.csv')
myco_tax<-read.csv('Processed_data/Bag_dat_myco_tax.csv')


# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from one or more groups of veg class
multipatt(Bag_Seq_wide %>% select(starts_with('ITSall')), # first argument is the community table, select only those columns
          Bag_Seq_wide$Fire.Severity) -> res.severity

multipatt(Bag_Seq_wide %>% select(starts_with('ITSall')), # first argument is the community table, select only those columns
          Bag_Seq_wide$Fire.Interval) -> res.interval

summary(res.interval)
summary(res.severity)

total_myco_reads<-Bag_Seq_wide%>% select(starts_with('ITSall'))%>% 
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
out.severity <- left_join(out.severity, Bag_Seq_wide %>% 
              select(Site,Transect,Location, Fire.Severity, starts_with('ITSall'))%>% 
              pivot_longer(cols = starts_with('ITSall'), names_to = 'OTU', values_to = 'count'))%>%
                group_by(OTU,Site,Transect,Location,Fire.Severity,count)%>%
                summarize(
                  count=sum(count))%>%
  mutate(OTU = fct_reorder(OTU, count, max), 
         Fire.Severity = as_factor(Fire.Severity),
         Location=as.factor(Location))%>%
  left_join(myco_tax)


out.interval <- left_join(out.interval, Bag_Seq_wide %>% 
              select(Site,Transect,Location,Fire.Interval,starts_with('ITSall'))%>% 
  pivot_longer(cols = starts_with('ITSall'), names_to = 'OTU', values_to = 'count'))%>%
  mutate(OTU = fct_reorder(OTU, count, max), 
         Interval = as_factor(Fire.Interval),
         Location=as.factor(Location))%>%
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
Bag_Seq_wide<-Bag_Seq_wide%>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

# extract the community table, save as a new object
mat_myco<- Bag_Seq_wide %>% select(starts_with("ITSall"))

# tax_gs_bag<-read.csv('Processed_data/taxa_w_genome_size_bag_data.csv')
# 
# myco_gs_otus <- tax_gs_bag$OTU
# 
# mat_ecm_gs<-mat_myco %>% 
#   dplyr::select(all_of(myco_gs_otus)) 

# run three permanovas, each with a different distance index / raw data input

adonis2(mat_ecm~ Fire.Severity+Fire.Interval + Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
         log10_biomass_day+myco_host_total 
          , data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)

table(Bag_Seq_wide$Fire.Interval)
table(Bag_Seq_wide$Fire.Severity)

# adonis2(mat_ecm_gs~ Fire.Severity+Fire.Interval + Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg +
#           log10_biomass_day+myco_host_total+gs_mean
#         ,data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)


cap.all <- capscale(mat_ecm~ Fire.Interval + Fire.Severity + Ortho_P_mg_kg+
                      Nitrate_mg_kg+ Ammonia_mg_kg+ perc_myco_host_freq + #myco_host_total
                      log10_biomass_day 
                    , data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)
anova(cap.all)
summary(cap.all)
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
cap.all
plot(cap.all)
proportions<-round(cap.all$CCA$eig/cap.all$tot.chi *100, 1) # proportion of variation associated with each axis

cap_test <- capscale(mat_ecm ~ 1 , data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)

ordistep(cap_test, formula(cap.all), direction='forward')



# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap.all, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site<- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')
scrs_biplot <- scrs %>% filter(score=='biplot')

interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_colors <- c("Hgih" = "blue", "Low" = "green")

library(ggrepel)

# first plot - site scores along with centroids for each group
cbind(Bag_Seq_wide, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2)) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_point(aes( colour= Fire.Interval,shape= Fire.Severity), size=8, stroke = 3)+ 
  geom_text(aes( label = Site), color= 'black', size=3)+
  #geom_text(data = scrs_cent, aes(label = label), size = 2) + 
  scale_shape_manual(values = c(19,1))+
  scale_colour_manual(values = interval_colors) +     # Custom colors for Interval
  labs( x=  paste0("CAP1 (", proportions[1], "%)"), y=  paste0("CAP2 (", proportions[2], "%)"))+
  geom_segment(data=scrs_biplot%>%
                 filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5),
               inherit.aes = FALSE,
               aes(x=0,y=0, xend=CAP1, yend=CAP2, group=label),
               arrow = arrow(type = "closed",length=unit(3,'mm')),
               color= 'black') +
  geom_text_repel(data=scrs_biplot%>%
                    filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5),
                  aes(x=CAP1, y=CAP2, label=label),
                  colour='black',size=6)+
  xlim(c(min(scrs_site[, 'CAP1']), max(scrs_site[, 'CAP1']))) + 
  ylim(c(min(scrs_site[, 'CAP2']), max(scrs_site[, 'CAP2']))) + 
  theme_bw() + 
  theme(legend.position='top')->p2

p2


# Perform NMDS

# Run NMDS (use distance = "robust.aitchison" if applicable)
set.seed(42)  # For reproducibility
nmds_model <- metaMDS(mat_ecm, distance = "robust.aitchison", k = 2, trymax = 100)

# Extract NMDS site scores
scrs <- as.data.frame(scores(nmds_model, display = "sites")) %>%
  mutate(Site = rownames(.), 
         Fire.Interval = Bag_Seq_wide$Fire.Interval, 
         Fire.Severity = Bag_Seq_wide$Fire.Severity)

# Extract biplot scores (environmental vectors)
biplot_scrs <- as.data.frame(scores(envfit(nmds_model, Bag_Seq_wide[, c("Ortho_P_mg_kg", "Nitrate_mg_kg", 
                                                                        "Ammonia_mg_kg", "perc_myco_host_freq", 
                                                                        "log10_biomass_day")]), display = "vectors")) %>%
  mutate(label = rownames(.))

# Define custom colors
interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_colors <- c("High" = "blue", "Low" = "green")

# Function to compute convex hulls
compute_hull <- function(data, group_var) {
  data %>%
    group_by({{ group_var }}) %>%
    slice(chull(NMDS1, NMDS2)) %>%
    ungroup()
}

# Compute convex hulls for Fire Interval and Fire Severity
hull_interval <- compute_hull(scrs, Fire.Interval)
hull_severity <- compute_hull(scrs, Fire.Severity)

# NMDS Plot with convex hulls
p_nmds <- ggplot(scrs, aes(x = NMDS1, y = NMDS2)) + 
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_polygon(data = hull_interval, aes(fill = Fire.Interval, group = Fire.Interval), 
               alpha = 0.3, color = "blue") +  # Convex hull for Fire Interval
  geom_polygon(data = hull_severity, aes(fill = Fire.Severity, group = Fire.Severity), 
               alpha = 0.3, color = "red") +  # Convex hull for Fire Severity
  geom_point(aes(color = Fire.Interval, shape = Fire.Severity), size = 6, stroke = 2) + 
  geom_text(aes(label = Site), color = 'black', size = 3) +
  scale_shape_manual(values = c(19, 1)) +
  scale_colour_manual(values = interval_colors) +  
  labs(x = "NMDS1", y = "NMDS2") +
  geom_segment(data = biplot_scrs,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2, group = label),
               arrow = arrow(type = "closed", length = unit(3, 'mm')),
               color = 'black') +
  geom_text_repel(data = biplot_scrs,
                  aes(x = NMDS1, y = NMDS2, label = label),
                  colour = 'black', size = 6) +
  theme_bw() + 
  theme(legend.position = 'top')

# Print the NMDS plot
p_nmds








# plot side-by-side using the patchwork package
library(patchwork)
p1 + p2


# top ten otus associated with the top of the ordination, presumablymost sensative samples
scrs_spp %>% 
  arrange(desc(CAP2)) %>% 
  head(10) %>% rownames_to_column(var= 'OTU')%>%
  left_join(tax)-> ind_otus
# taxonomic information for those otus
tax %>% 
  filter(OTU %in% ind_otus$label)


# still tidying - plotting turnover in space
pco1 <- capscale(mat_ecm ~ 1, data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)
scrs_site <- scores(pco1, display='sites') # TIDY
cbind(Bag_Seq_wide, scrs_site) %>% 
  ggplot(aes(x=Longitude, y=Latitude)) + 
  geom_point()


# including spatial patterns in analyses of community composition
# one way to do this is using principle coordinates of neighbour matrices
xy <- dist(Bag_Seq_wide[, c('Longitude', 'Latitude')]) # create distance matrix of longs and lats
xy.pcnm <- pcnm(xy)$vectors %>% as.data.frame() # export the result into a dataframe

# combine the PCNM results with the rest of the data  
temp <- cbind(Bag_Seq_wide, xy.pcnm)

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
 ordistep(cap.0, formula(cap.sp), direction='forward')


# should we include all climate variables in our analysis
# check for variance inflation - values higher than ~ 10 are unlikely to explain unique variation
cap.cl <- capscale(mat_ecm ~ Fire.Interval + Fire.Severity + Ortho_P_mg_kg+
                     Nitrate_mg_kg+ Ammonia_mg_kg+ perc_myco_host_freq + #myco_host_total
                     log10_biomass_day , data=Bag_Seq_wide)

vif.cca(cap.cl)


# test unique variation explained by partitions using partial CAPs
# 1 - associated with Fire
anova(capscale(mat_ecm ~ Fire.Interval+Fire.Severity +
                 Condition(Longitude + Latitude+ Ortho_P_mg_kg+
                             Nitrate_mg_kg+ Ammonia_mg_kg+ perc_myco_host_freq + #myco_host_total
                             log10_biomass_day  ), data=temp, 
               distance='robust.aitchison'))
# 2 - associated with soil nutrients
anova(capscale(mat_ecm ~ Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg+
                 Condition(Fire.Interval + Fire.Severity+  perc_myco_host_freq + #myco_host_total
                             log10_biomass_day ), data=temp, 
               distance='robust.aitchison'))
# 3 - associated with space
anova(capscale(mat_ecm ~ Longitude + Latitude +
                 Condition(Fire.Interval+Fire.Severity + Ortho_P_mg_kg+
                  Nitrate_mg_kg+ Ammonia_mg_kg+ perc_myco_host_freq), data=temp, 
               distance='robust.aitchison'))
# variation partitioning - here to three groups of variables
vp <- varpart(vegdist(mat_ecm, distance='robust.aitchison'), 
              ~ Fire.Interval+ Fire.Severity, # Influence of Fire Regime = X1
              ~Ortho_P_mg_kg+Nitrate_mg_kg+ Ammonia_mg_kg, #nutes
              ~ perc_myco_host_freq , # myco host
              #  ~ PCNM2, # spatial = X3 (I selected this variable from the ordistep() lines 334:337)
              ~ Site #4 effect of site
              , data=temp)
vp
plot(vp, Xnames=c('Regime','Nutrients','Myco_host','Site'))




