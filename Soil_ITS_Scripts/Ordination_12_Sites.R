
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
Bag_Site<-read.csv('Processed_data/All_Bag_Site_Info.csv')

Nutrients_Transects<-read_excel('Processed_data/Nutrients_Transect_level.xlsx', col_types = 'numeric')
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


Blast_ID<-Blast_ID%>%
  dplyr::mutate(Regime = paste(Interval, Severity, sep = "_"),
                Site=as.numeric(Site),
                transect=as.numeric(transect))%>%
  rename(Transect=transect)%>%
  #add the pairs of thee site into the df
  left_join(Bag_Site %>% select(Site,Transect, Site_Pair)%>% unique(), by = c("Site","Transect"))%>%
  left_join(VEG_COVER_Transects%>%
              mutate(Site=as.numeric(Site), Transect=as.numeric(Transect)))%>%
  left_join(Nutrients_Transects)%>%
  select(-`Fire 3`,-`Interval (yrs)...16`,-`FESM severity category`)



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

tax <- tax %>%
  mutate(Species = ifelse(str_detect(Species, '_sp'), NA, Species),
         across(everything(), ~ str_remove(., "_fam_Incertae_sedis|_gen_Incertae_sedis|_ord_Incertae_sedis|_cls_Incertae_sedis|_phy_Incertae_sedis")))


library(readxl)
Fun_Traits<-read_excel("Processed_data/Polme_etal_2021_FungalTraits.xlsx")
# add exploration types
tax<-tax%>%
  left_join(Fun_Traits%>%rename(Genus=GENUS))%>%
  select(SH_ID:guild,Ectomycorrhiza_exploration_type_template,Ectomycorrhiza_lineage_template)%>%
  rename(exploration_type=Ectomycorrhiza_exploration_type_template,
         Ecm_lineage=Ectomycorrhiza_lineage_template)%>%
  mutate(across(c(exploration_type,Ecm_lineage), ~replace_na(.x, "unknown")))# Replace NA with "Unknown" 






# transpose community table, to sample-taxon, for vegan functions and remove taxonomy col
otu%>%
  select(-last_col())%>%
  remove_rownames()%>%
  column_to_rownames("SH_ID")%>%
  t()-> mat


# assess variation in sampling effort, plotting sample effort curves
temp <- rarecurve(mat, step=1000, tidy=TRUE)
Blast_ID%>%
  left_join( temp %>% rename( sample_ID=Site))%>% 
  #filter(Site %in% c( 50:63))%>%
  ggplot(aes(x=Sample, y=Species, colour=as.factor(Transect), group=sample_ID)) + 
  geom_line() 

# rarefy the community matrix, using an arbitrary cut-off
mat<-rrarefy(mat, 7000)
matr <- temp[rowSums(mat)==7000, ]
dim(matr)


# prepare our data for betadiversity analyses:
# 1 - focus analysis on ectomycorrhizal fungal taxa (create vector of OTU ids)
# 2 - join the community table with our metadata table (need to also modify the community table so that it can be joined)
ecm_otus <- filter(tax, guild=='ectomycorrhizal')$SH_ID
dat_ecm <- left_join(Blast_ID,  mat %>% 
            as.data.frame() %>% 
            dplyr::select(ecm_otus) %>% # just ecto OTUs
            rownames_to_column('sample_ID'))

#this is now only looking at the 12 sites I have selected
dat_ecm_12_site<-dat_ecm%>% filter(!is.na(Site_Pair))

dat_ecm_12_site %>%
  count(Severity, Interval)

# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from one or more groups of veg class
multipatt(dat_ecm_12_site %>% select(starts_with('SH')), # first argument is the community table, select only those columns
          dat_ecm_12_site$Interval) -> res
summary(res)

# to visualise differences in taxonomic composition, using output from multipatt()
# first prepare the data in the res object so that it can be joined with taxonomic info
# output is only those otus significant for one or more pairs
out <- res[['sign']] %>% #what does this do?
  filter(p.value <= 0.05) %>% 
  rownames_to_column('SH_ID') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(site_code = gsub('^s.', '', group))
# then join with the taxonomy table, then the relevant community data, 
# and reorder the otu levels by decreasing abundance
out.Interval <- left_join(out, tax) %>% 
  left_join(dat_ecm_12_site %>% 
              select(Site,Transect, sample,sample_ID, Interval, starts_with('SH')) %>% 
              pivot_longer(cols=starts_with('SH'), names_to='SH_ID', 
                           values_to='count')) %>% 
  mutate(OTU_ID = fct_reorder(SH_ID, count, max), 
         Interval = as_factor(Interval))



#calculate relative abundance
p<-dat_ecm_12_site%>% select(ends_with('.09FU'))%>% 
  pivot_longer(cols=ends_with('.09FU'), names_to='SH_ID', 
               values_to='count')%>%left_join(tax)%>%
  summarise(sum(count))%>%
  pull()

p

# finally produce the barplot
p1<-out.Interval %>% 
  ggplot(aes(x=Interval, y=(count/p)*100, fill=exploration_type, text=SH_ID)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
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
p1

# use this next lines to interactively get OTU IDs
plotly::ggplotly(p1)
#######################


###next analysis#######

# next analysis - permanova
# extract the community table, save as a new object
mat_ecm <- dat_ecm_12_site %>% select(ends_with('.09FU'))

dat_ecm_12_site<-dat_ecm_12_site%>%
  mutate(
    Bray.P = ifelse(is.na(Bray.P), mean(Bray.P, na.rm = TRUE), Bray.P),
    Total.P = ifelse(is.na(Total.P), mean(Total.P, na.rm = TRUE), Total.P),
    NO3 = ifelse(is.na(NO3), mean(NO3, na.rm = TRUE), NO3) )



#remove zeros
 rows_zero<-which(apply(mat_ecm, 1, function(row) all(row == 0)))
# 
df_rows_zero<-dat_ecm_12_site[rows_zero,]
# 
 mat_ecm_w<-mat_ecm[-rows_zero,]
# 
dat_ecm_12_site_w<-dat_ecm_12_site[-rows_zero,]

# run three permanovas, each with a different distance index / raw data input

adonis2(mat_ecm_w ~ Severity+Interval+ Veg_Class +
          Tree.Basal.Area_m2 + Herb.Cover_0.50cm_perc +
          NH4 + NO3 + Total.P, data=dat_ecm_12_site_w, distance='robust.aitchison', add=TRUE)

table(dat_ecm_12_site_w$Interval)



Nute_Veg<-dat_ecm_12_site%>%
  select(Tree.Basal.Area_m2:Nitrogen,Fire.Interval,Fire.Severity)

rows_zero_na <- which(apply(Nute_Veg, 1, function(row) all(row == 0) | any(is.na(row))))

dat_ecm_12<-dat_ecm_12_site%>%
  select(Tree.Basal.Area_m2:SH1205826.09FU)
 colnames(dat_ecm_12_site)

cap.all <- capscale(dat_ecm_12~ Interval+Severity +
                      NO3 + Total.P  + #NH4
                      Shrub.Cover_50.200cm_perc + Tree.Basal.Area_m2
                    , data=dat_ecm_12_site, distance='robust.aitchison', add=TRUE)
anova(cap.all, by = "margin")
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
#write_xlsx(Cap1_aov,'Processed_Data/CAP_12_Sites_AOV.xlsx') 
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

# first plot - site scores along with centroids for each group
cbind(dat_ecm_12_site, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2)) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_point(aes( colour= Interval ,shape = Severity), size=6)+ 
  geom_text(aes( label = label), color= 'black', size=3)+
  #geom_text(data = scrs_cent, aes(label = label), size = 2) + 
  scale_colour_manual(values = interval_colors) +     # Custom colors for Interval
  labs( x=  paste0("CAP1 (", proportions[1], "%)"), y=  paste0("CAP2 (", proportions[2], "%)"))+
  # geom_segment(data=scrs_spp%>%
  #                filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5),
  #              inherit.aes = FALSE,
  #              aes(x=0,y=0, xend=CAP1, yend=CAP2, group=label),
  #              arrow = arrow(type = "closed",length=unit(3,'mm')),
  #              color= 'black') +
  # geom_text_repel(data=scrs_spp%>%
  #                   filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5)%>%
  #                   rename(SH_ID=label)%>%
  #                   left_join(tax),
  #                 inherit.aes = FALSE,
  #                 aes(x=CAP1, y=CAP2, label=Genus),
  #                 colour='black',size=6)+
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
  head(10) -> ind_otus
# taxonomic information for those otus
tax %>% 
  filter(SH_ID %in% ind_otus$label)


# still tidying - plotting turnover in space
pco1 <- capscale(dat_ecm_12 ~ 1, data=dat_ecm_12_site, distance='robust.aitchison', add=TRUE)
scrs_site <- scores(pco1, display='sites') # TIDY
cbind(dat_ecm_12_site, scrs_site) %>% 
  ggplot(aes(x=Longitude, y=Latitude)) + 
  geom_point()


# including spatial patterns in analyses of community composition
# one way to do this is using principle coordinates of neighbour matrices
xy <- dist(dat_ecm_12_site[, c('Longitude', 'Latitude')]) # create distance matrix of longs and lats
xy.pcnm <- pcnm(xy)$vectors %>% as.data.frame() # export the result into a dataframe

# combine the PCNM results with the rest of the data  
temp <- cbind(dat_ecm_12_site, xy.pcnm)

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
cap.cl <- capscale(dat_ecm_12 ~ Interval+Severity +
                     NO3 + Total.P  + #NH4
                     Shrub.Cover_50.200cm_perc + Tree.Basal.Area_m2, data=dat_ecm_12_site)

vif.cca(cap.cl)

# variation partitioning - here to three groups of variables
vp <- varpart(vegdist(dat_ecm_12, distance='robust.aitchison'), 
              ~Interval+Severity, # type of environment sample was collected from = X1
              #~elev , # climate = X2
              ~ PCNM2 + PCNM9  # spatial = X3
            , data=temp)
vp
plot(vp, Xnames=c('Interval', 'Climate', 'Space'))

# test unique variation explained by partitions using partial CAPs
# 1 - associated with site_code
anova(capscale(mat_ecm ~ Interval +
                 Condition(Annual_Temp + elev +
                              PCNM2 + PCNM9 ), data=temp, 
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





