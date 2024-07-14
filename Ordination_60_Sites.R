
#LOAD LIBRARARIES
library(forcats)
library(tidyr)
library(dplyr)
library(vegan)
library(indicspecies)
library(tibble) 
library(readxl)
library(ggplot2)
library(ggrepel)


#Enviromental data
precip<-read_excel('Processed_data/Site_Precip.xlsx')

Precip_Site<-precip%>%
  select(Site,chirps,MonthYear)%>%
  group_by(MonthYear,Site)%>%
  summarise(Rain=mean(chirps), .groups = 'keep')%>%
  distinct()%>%
  ungroup()%>%
  select(-MonthYear)%>%
  group_by(Site)%>%
  summarise(Rain=mean(Rain), .groups = 'keep')%>%
  mutate(Site=gsub('ABS0|ABS00','',Site))


Meta_Site <- read_excel("Processed_data/Site_Precip_Temp_Elv.xlsx")

Meta_Site<-Meta_Site%>%
  select(wc2.1_30s_elev,Annual_Temp,site)%>%
  rename(Site=site,
         elev= wc2.1_30s_elev)%>%
  mutate(Site=gsub('ABS0|ABS00','',Site))%>%
 left_join(Precip_Site)


#All meta data from 12 sites with bags collected
Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')

Nutrients_Transects<-read_excel('Processed_data/Nutrients_Transect_level.xlsx')
VEG_COVER_Transects <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                  sheet = "Transect.Level_Data")
VEG_COVER_Transects$Site= sub(c('ABS00|ABS0'),'',VEG_COVER_Transects$Site)
VEG_COVER_Transects$Transect= sub(c('T'),'',VEG_COVER_Transects$Transect)


#import Blast ID df and clean data
Blast_ID<-read_excel('Raw_data/Updated_Data/ABS.MER.fielddata.Feb.2023_Site.Info_AF.xlsx')
Blast_ID$Site= sub(c('ABS00|ABS0'),'',Blast_ID$Site)
Blast_ID$sample_ID<-as.character(Blast_ID$sample_ID)
Blast_ID$transect<-as.character(Blast_ID$transect)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
names(Blast_ID)[9]<-'Veg_Class_Abv'

Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')


Blast_ID<-Blast_ID%>%
  dplyr::mutate(Regime = paste(Interval, Severity, sep = "_"))%>%
  rename(Transect=transect)%>%
  #add the pairs of thee site into the df
  left_join(Bag_Site %>% select(Site,Transect, Site_Pair)%>% unique(), by = c("Site","Transect"))%>%
  left_join(VEG_COVER_Transects)%>%
  select(-`Fire 3`,-`Interval (yrs)...16`,-`FESM severity category`)%>%
  left_join(Meta_Site)


#funguild output
AM <- read_excel("Funguild/AM.xlsx")
Ecto <- read_excel("Funguild/Ecto.xlsx")
Path<- read_excel("Funguild/Path.xlsx")
Sap<- read_excel("Funguild/Sap.xlsx")


#taxonomy table
otu<-read.csv('Raw_data/Updated_Data/demultiplexed.cleaned.combined.cf.fasta.blast.i97.a95.csv')
#remove first three rows that do not have taxonomy ids and 2 rows with totals and sample counts 
otu<-otu[-c(1:3),-c(2:3) ]

tax <- otu %>%
  separate(taxonomy, into = c("Percentage", "Species_Name", "Accession", "SH_ID", "Reference", "Taxonomy_Details"), sep = "\\|", extra = "drop") %>%
  separate(Taxonomy_Details, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  select(SH_ID:Species)%>%
  mutate(across(Kingdom:Species, ~ gsub('^[a-z]__', '', .)),
         guild = case_when(Genus %in% Ecto$Genus ~ 'ectomycorrhizal', 
                           Genus %in% AM$Genus ~ 'arbuscular_mycorrhizal', 
                           Genus %in% Path$Genus ~ 'pathogenic', 
                           Genus %in% Sap$Genus ~ 'saprotrophic'))
# any kingdom-to-genus unclassified, change to NA
tax[tax %in% grep('unclassified', tax, value=T)] <- NA
# any species unclassified, change to NA
tax[tax %in% grep('_sp$', tax, value=T)] <- NA

#change first col name to something easy
names(otu)[1]<-'SH_ID'

# transpose community table, to sample-taxon, for vegan functions and remove taxonomy col
mat<-otu%>%
  select(-last_col())%>%
  remove_rownames()%>%
  column_to_rownames("SH_ID")%>%
  t()

# assess variation in sampling effort, plotting sample effort curves
temp <- rarecurve(mat, step=1000, tidy=TRUE)
Blast_ID%>%
  left_join( temp %>% rename( sample_ID=Site))%>% 
  #filter(Site %in% c( 50:63))%>%
  ggplot(aes(x=Sample, y=Species, colour=as.factor(Transect), group=sample_ID)) + 
  geom_line() + 
  facet_wrap(~as.numeric(Site))



# check variation in sample effort, looking for break points WHAT DOES THIS MEAN
hist(log10(rowSums(mat)))
sort(rowSums(mat))[1:63]


# rarefy the community matrix, using an arbitrary cut-off
#Im not going to do this here though because it isnt needed
# temp<-rrarefy(mat, 5000)
# matr <- temp[rowSums(temp)==5000, ]
# dim(matr)


# prepare our data for betadiversity analyses:
# 1 - focus analysis on ectomycorrhizal fungal taxa (create vector of OTU ids)
# 2 - join the community table with our metadata table (need to also modify the community table so that it can be joined)
ecm_otus <- filter(tax, guild=='ectomycorrhizal')$SH_ID
dat_ecm <- left_join(Blast_ID,  mat %>% 
                       as.data.frame() %>% 
                       dplyr::select(ecm_otus) %>% # just ecto OTUs
                       rownames_to_column('sample_ID'))


# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from fire interval
res<-multipatt(dat_ecm%>% select(starts_with('SH')), # first argument is the community table, select only those columns
               dat_ecm$Regime) 
summary(res)

# to visualise differences in taxonomic composition, using output from multipatt()
# first prepare the data in the res object so that it can be joined with taxonomic info
# output is only those otus significant for one or more pairs
out <- res[['sign']] %>% #what does this do?
  filter(p.value <= 0.05) %>% 
  rownames_to_column('SH_ID') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(Regime = gsub('^s.', '', group))
# then join with the taxonomy table, then the relevant community data, 
# and reorder the otu levels by decreasing abundance
out <- left_join(out, tax) %>% 
  left_join(dat_ecm %>% 
              select(Site,Transect, sample,sample_ID, Regime, ends_with('.09FU')) %>% 
              pivot_longer(cols=ends_with('.09FU'), names_to='SH_ID', 
                           values_to='count')) %>% 
  mutate(SH_ID = fct_reorder(SH_ID, count, max), 
         Regime = as_factor(Regime))

library(RColorBrewer)


# Generate a custom color palette by combining multiple RColorBrewer palettes
custom_palette <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"), brewer.pal(9, "Set1"))

# Ensure you have enough unique colors
custom_palette <- unique(custom_palette)

# finally produce the barplot
out %>% 
  ggplot(aes(x=Regime, y=count, fill=Genus, text=SH_ID)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat='identity', position=position_fill()) + 
  scale_x_discrete(drop=FALSE) + 
  scale_fill_manual(values = custom_palette) +  #palette.pals()
  scale_y_continuous(labels = scales::percent) + 
  labs(y='Percentage') + 
  theme_bw() -> p1
p1

# use this next lines to interactively get OTU IDs
plotly::ggplotly(p1)
#######################


####next analysis#######

# next analysis - permanova
# extract the community table, save as a new object
mat_ecm <- dat_ecm %>% select(ends_with('.09FU'))


# 
rows_zero_na <- which(apply(mat_ecm, 1, function(row) all(row == 0 | is.na(row))))
# 
zero_value_rows<-dat_ecm[rows_zero_na,]
# 
mat_ecm<-mat_ecm[-rows_zero_na,]
# 
dat_ecm<-dat_ecm[-rows_zero_na,]



# Identify rows that contain any NA values in the specified columns
rows_with_na <- which(apply(dat_ecm, 1, function(row) any(is.na(row))))

# Print the rows with NA values
na_rows <- dat_ecm[rows_with_na, ]


#this is to look at all 60 sites
#above is to only look at 30 sites
#mat_ecm <- dat_ecm %>% select(ends_with('.09FU'))



adonis2(mat_ecm ~ Severity*Interval, data=dat_ecm, distance='robust.aitchison', add=TRUE)

table(dat_ecm$Interval)


# a constrained analysis of principal coordinates using a different distance index - result is quite good
cap1 <- capscale(mat_ecm ~ Interval*Severity +
                   Shrub.Cover_50.200cm_perc + Tree.Basal.Area_m2 +
                   Condition (Site)
                 , data=dat_ecm, distance='robust.aitchison', add=TRUE)
Cap1_aov<-as.data.frame( anova(cap1, by = "margin"))%>%
  rownames_to_column()
#write_xlsx(Cap1_aov,'Processed_Data/CAP_60_Sites_AOV.xlsx')
plot(cap1)
cap1 # summary of inertia
proportions<-round(cap1$CCA$eig/cap1$tot.chi *100, 1) # proportion of variation associated with each axis
anova(cap1) # statistical significance of the constrait

cap_test <- capscale(mat_ecm ~ 1 , data=dat_ecm, distance='robust.aitchison', add=TRUE)

ordistep(cap_test, formula(cap_test), direction='forward')



# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap1, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site <- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')
scrs_biplot <- scrs %>% filter(score=='biplot')


interval_colors <- c("Long" = "darkred", "Short" = "orange")

# first plot - site scores along with centroids for each group
cbind(dat_ecm, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=MDS1)) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_point(aes( colour= Interval ,shape = Severity), size=6)+ 
  geom_text(aes( label = label), color= 'black', size=3)+
  #geom_text(data = scrs_cent, aes(label = label), size = 2) + 
  scale_colour_manual(values = interval_colors) +     # Custom colors for Interval
  #labs( x=  paste0("CAP1 (", proportions[1], "%)"), y=  paste0("MDS1 (", proportions[2], "%)"))+
  geom_segment(data=scrs_spp%>%
                 filter(abs(CAP1) > 0.5 | abs(MDS1) > 0.5),
               inherit.aes = FALSE,
               aes(x=0,y=0, xend=CAP1, yend=MDS1, group=label),
               arrow = arrow(type = "closed",length=unit(3,'mm')),
               color= 'black') +
  geom_text_repel(data=scrs_spp%>%
                    filter(abs(CAP1) > 0.5 | abs(MDS1) > 0.5)%>%
                    rename(SH_ID=label)%>%
                    left_join(tax),
                  inherit.aes = FALSE,
                  aes(x=CAP1, y=MDS1, label=Genus),
                  colour='black',size=7)+
  xlim(c(min(scrs_site[, 'CAP1']), max(scrs_site[, 'CAP1']))) + 
  ylim(c(min(scrs_site[, 'MDS1']), max(scrs_site[, 'MDS1']))) + 
  theme_bw() + 
  theme(legend.position='top')->p2

p2

# plot side-by-side using the patchwork package
library(patchwork)
p1 + p2


# top ten otus associated with the top of the ordination, presumably '???' samples
scrs_spp %>% 
  arrange(desc(CAP2)) %>% 
  head(10) -> inv_otus
# taxonomic information for those otus
tax %>% 
  filter(SH_ID %in% inv_otus$label)


# still tidying - plotting turnover in space
pco1 <- capscale(mat_ecm ~ 1, data=dat_ecm, distance='robust.aitchison', add=TRUE)
scrs_site <- scores(pco1, display='sites') # TIDY
cbind(dat_ecm, scrs_site) %>% 
  ggplot(aes(x=Longitude, y=Latitude, colour=MDS1)) + 
  geom_point()


# including spatial patterns in analyses of community composition
# one way to do this is using principle coordinates of neighbour matrices
xy <- dist(dat_ecm[, c('Longitude', 'Latitude')]) # create distance matrix of longs and lats
xy.pcnm <- pcnm(xy)$vectors %>% as.data.frame() # export the result into a dataframe

# combine the PCNM results with the rest of the data  
temp <- cbind(dat_ecm, xy.pcnm)

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
cap.cl <- capscale(mat_ecm ~ Annual_Temp + elev, data=dat_ecm)
vif.cca(cap.cl)

# variation partitioning - here to three groups of variables
vp <- varpart(vegdist(mat_ecm, distance='robust.aitchison'), 
              ~Interval, # type of environment sample was collected from = X1
              ~ elev , # climate = X2
              ~ PCNM2 + PCNM9  # spatial = X3
              , data=temp)
vp
plot(vp, Xnames=c('Interval', 'Climate', 'Space'))

# test unique variation explained by partitions using partial CAPs
# 1 - associated with site_code
anova(capscale(mat_ecm ~ Interval +
                 Condition(elev +
                             PCNM2 + PCNM9 ), data=temp, 
               distance='robust.aitchison'))
# 2 - associated with climate
anova(capscale(mat_ecm ~ elev +
                 Condition(Interval + 
                             PCNM2 + PCNM9), data=temp, 
               distance='robust.aitchison'))
# 3 - associated with space
anova(capscale(mat_ecm ~ PCNM2 + PCNM9 +
                 Condition( elev + Interval), data=temp, 
               distance='robust.aitchison'))




