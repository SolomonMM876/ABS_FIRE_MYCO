
#LOAD LIBRARARIES
library(forcats)
library(tidyr)
library(dplyr)
library(vegan)
library(indicspecies)
library(tibble) 
library(readxl)
library(ggplot2)




#All meta data from 12 sites with bags collected
Bag_data<-read.csv('Processed_data/All_Bag_data.csv')[,-1]
  
Seq_data<-read.csv('Processed_data/cleaned_seq_dat.csv')%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))%>%
  filter(guild %in%c("Ectomycorrhizal","Arbuscular Mycorrhizal","Orchid Mycorrhizal","Ericoid Mycorrhizal"))%>%
  select(-sample)


tax<-Seq_data%>%
  select(OTU,kingdom:species, guild)%>%
  distinct()

#format required for some vegan functions
wide_seq<-Seq_data  %>%
    pivot_wider(
    names_from = OTU,         # Column containing OTU names that will become new columns
    values_from = count,      # Column containing values that will fill the new columns
    id_cols = Tube_ID,  # Column(s) to keep as identifier
    values_fill = 0
    )

Bag_Seq_wide<-left_join(Bag_data%>%
           select(Tube_ID,Fire.Interval,Fire.Severity,Site, Transect, Location,#Site Data
                  Ortho_P_mg_kg,Nitrate_mg_kg, Ammonia_mg_kg,#Nutrients from resin data
                  log10_biomass_day, #biomass data
                  Longitude,Latitude), #meta data
                        wide_seq)%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))

#check to make sure there are no NAs
NA_rows <- Bag_Seq_wide %>%
  mutate(across(everything(), ~ ifelse(is.na(.), ., NA))) %>%
  filter(if_any(everything(), ~ !is.na(.)))


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
               values_to='count')%>%left_join(tax)%>%
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
  left_join(tax)


out.interval <- left_join(out.interval, Bag_Seq_wide %>% 
              select(Site,Transect,Location,Fire.Interval,starts_with('ITSall'))%>% 
  pivot_longer(cols = starts_with('ITSall'), names_to = 'OTU', values_to = 'count'))%>%
  mutate(OTU = fct_reorder(OTU, count, max), 
         Interval = as_factor(Fire.Interval),
         Location=as.factor(Location))%>%
  left_join(tax)

# then join with the taxonomy table
tax.interval <- left_join(out.interval, tax)


# finally produce the barplot
out.severity.filter<-out.severity%>%
  filter(count>0)%>%
  mutate( rel_abun= (count/total_myco_reads)*100)

out.severity.filter %>% 
  ggplot(aes(x=Fire.Severity, y=rel_abun, fill=family, color=species, text=OTU)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat='identity', position=position_fill(), color='black') + 
  scale_x_discrete(drop=FALSE) + 
  scale_fill_brewer(palette='Set3') + 
  # scale_y_continuous(labels = scales::percent) + 
  labs(y='Percentage') + 
  theme_bw()+
  facet_wrap(~ Site)-> p1
p1


out.interval<-out.interval %>% 
  filter(count>0)%>%
  mutate( rel_abun= (count/total_myco_reads)*100)
  



# use this next lines to interactively get OTU IDs
plotly::ggplotly(p1)
#######################


###next analysis#######

# next analysis - permanova

# first remove samples that have no mycorrhizal reads because they cause errors below
Bag_Seq_wide<-Bag_Seq_wide%>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

# extract the community table, save as a new object
mat_ecm<- Bag_Seq_wide %>% select(starts_with("ITSall"))



# run three permanovas, each with a different distance index / raw data input

adonis2(mat_ecm~ Fire.Severity+Fire.Interval + Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg,
        data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)

table(Bag_Seq_wide$Fire.Interval)
table(Bag_Seq_wide$Fire.Severity)




cap.all <- capscale(mat_ecm~ Fire.Interval + Fire.Severity + Ortho_P_mg_kg+
                      Nitrate_mg_kg+ Ammonia_mg_kg+
                      log10_biomass_day+
                      Condition (Site)
                    , data=Bag_Seq_wide, distance='robust.aitchison', add=TRUE)
anova(cap.all, by = "margin")
summary(cap.all)
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
cap.all
plot(cap.all)
proportions<-round(cap.all$CCA$eig/cap.all$tot.chi *100, 1) # proportion of variation associated with each axis


# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap.all, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site<- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')

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
  geom_segment(data=scrs_spp%>%
                 filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5),
               inherit.aes = FALSE,
               aes(x=0,y=0, xend=CAP1, yend=CAP2, group=label),
               arrow = arrow(type = "closed",length=unit(3,'mm')),
               color= 'black') +
  geom_text_repel(data=scrs_spp%>%
                    filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5)%>%
                    rename(OTU=label)%>%
                    left_join(tax),
                  inherit.aes = FALSE,
                  aes(x=CAP1, y=CAP2, label=genus),
                  colour='black',size=6)+
  xlim(c(min(scrs_site[, 'CAP1']), max(scrs_site[, 'CAP1']))) + 
  ylim(c(min(scrs_site[, 'CAP2']), max(scrs_site[, 'CAP2']))) + 
  theme_bw() + 
  theme(legend.position='top')->p2

p2




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
cap.cl <- capscale(mat_ecm ~ Ortho_P_mg_kg+
                     Nitrate_mg_kg+ Ammonia_mg_kg, data=Bag_Seq_wide)

vif.cca(cap.cl)


# test unique variation explained by partitions using partial CAPs
# 1 - associated with Fire
anova(capscale(mat_ecm ~ Fire.Interval+Fire.Severity +
                 Condition(Longitude + Latitude ), data=temp, 
               distance='robust.aitchison'))
# 2 - associated with climate
anova(capscale(mat_ecm ~ Annual_Temp + elev +
                 Condition(Interval + 
                             PCNM2 + PCNM9), data=temp, 
               distance='robust.aitchison'))
# 3 - associated with space
anova(capscale(mat_ecm ~ PCNM2 + PCNM9 +
                 Condition(Annual_Temp + elev + Interval), data=temp, 
               distance='robust.aitchison'))





