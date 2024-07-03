
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
Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')

Bag_Site<-left_join(Bag_Site, Bag_Site%>%
   filter(!is.na(Pair))%>%
   select(Site,Pair)%>%
  unique()%>%
  group_by(Pair) %>%
  summarize(Site_Pairs = paste(Site, collapse = "_")))


#import Blast ID df and clean data
Blast_ID<-read_excel('Raw_data/Updated_Data/ABS.MER.fielddata.Feb.2023_Site.Info_AF.xlsx')
Blast_ID$Site= sub(c('ABS00|ABS0'),'',Blast_ID$Site)
Blast_ID$sample_ID<-as.character(Blast_ID$sample_ID)
Blast_ID$transect<-as.character(Blast_ID$transect)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
names(Blast_ID)[9]<-'Veg_Class_Abv'

Blast_ID<-Blast_ID%>%
  dplyr::mutate(Regime = paste(Interval, Severity, sep = "_"))%>%
  rename(Transect=transect)%>%
  #add the pairs of thee site into the df
  left_join(Bag_Site %>% select(Site,Transect, Site_Pairs)%>% unique(), by = c("Site","Transect"))


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

#this is now only looking at the 12 sites I have selected
dat_ecm_12_site<-dat_ecm%>% filter(!is.na(Site_Pairs))


# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from one or more groups of veg class
multipatt(dat_ecm_12_site %>% select(starts_with('SH')), # first argument is the community table, select only those columns
          dat_ecm_12_site$Site_Pairs) -> res
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
out <- left_join(out, tax) %>% 
  left_join(dat_ecm_12_site %>% 
              select(Site,Transect, sample,sample_ID, Site_Pairs, starts_with('SH')) %>% 
              pivot_longer(cols=starts_with('SH'), names_to='SH_ID', 
                           values_to='count')) %>% # , relationship = "many-to-many" WARNing message here
  mutate(OTU_ID = fct_reorder(SH_ID, count, max), 
         Site_Pairs = as_factor(Site_Pairs))

# finally produce the barplot
out %>% 
  ggplot(aes(x=Site_Pairs, y=count, fill=Order, text=SH_ID)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat='identity', position=position_fill()) + 
  scale_x_discrete(drop=FALSE) + 
  scale_fill_brewer(palette='Set3') + 
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
mat_ecm <- dat_ecm %>% select(starts_with('SH'))

rows_all_zero<-which(apply(mat_ecm, 1, function(row) all(row == 0)))

zero_value_rows<-mat_ecm[rows_all_zero,]

mat_ecm<-mat_ecm[-rows_all_zero,]

dat_ecm<-dat_ecm[-rows_all_zero,]

# run three permanovas, each with a different distance index / raw data input
adonis2(mat_ecm ~ Site, data=dat_ecm, distance='bray', add=TRUE)
adonis2(mat_ecm ~ sample, data=dat_ecm, distance='bray', binary=TRUE, add=TRUE)
adonis2(mat_ecm ~ sample, data=dat_ecm, distance='robust.aitchison', add=TRUE)

table(dat_ecm$Site)
# for more distance/similarity indices, look at:
?dist
?vegdist
?labdsv::dsvdis

# for ways to standardise/transform the data prior to analysis, look at:
?decostand

# a redundancy analysis (constrained PCA) - ordination result is poor, don't interpret this
rda1 <- rda(mat_ecm ~ dat_ecm$sample, data=dat_ecm)
plot(rda1)

# a constrained analysis of principal coordinates using bray-curtis distances - ordination result is better, still some skew on the top of the y-axis
cap1 <- capscale(mat_ecm ~ dat_ecm$sample, data=dat_ecm, distance='bray', add=TRUE)
plot(cap1)

# try the same but with presence-absence instead of raw counts - still some skew but better
cap1 <- capscale(decostand(mat_ecm, method='pa') ~ 
                   dat_ecm$sample, data=dat_ecm, distance='bray', add=TRUE)
plot(cap1)

# try a different distance index - result is quite good
cap1 <- capscale(mat_ecm ~ dat_ecm$sample, data=dat_ecm, distance='robust.aitchison', add=TRUE)
plot(cap1)
cap1 # summary of inertia
cap1$CCA$eig/cap1$tot.chi # proportion of variation associated with each axis
anova(cap1) # statistical significance of the constrait

# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap1, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site <- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')

# first plot - site scores along with centroids for each group
cbind(dat_ecm, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2, colour=site_code)) + 
  geom_point(size=0.5, alpha=0.5) + 
  geom_point(data=scrs_cent %>% 
               rename(site_code=label) %>% 
               mutate(site_code=gsub('site_code', '', site_code)), size=2) + 
  theme_bw() + 
  theme(legend.position='top') -> p1
# second plot - species scores for those loaded heavily along at least one axis
# this one still needs some work
ggplot(scrs_spp %>% 
         filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5), aes(x=CAP1, y=CAP2, label=label)) + 
  geom_text() + 
  xlim(c(min(scrs_site[, 'CAP1']), max(scrs_site[, 'CAP1']))) + 
  ylim(c(min(scrs_site[, 'CAP2']), max(scrs_site[, 'CAP2']))) + 
  theme_bw() -> p2
# plot side-by-side using the patchwork package
library(patchwork)
p1 + p2


# top ten otus associated with the top of the ordination, presumably 'INV' samples
scrs_spp %>% 
  arrange(desc(CAP2)) %>% 
  head(10) -> inv_otus
# taxonomic information for those otus
tax %>% 
  filter(OTU_ID %in% inv_otus$label)


# still tidying - plotting turnover in space
pco1 <- capscale(mat_ecm1 ~ 1, data=dat_ecm, distance='bray', add=TRUE)
scrs_site <- scores(pco1, display='sites') # TIDY
cbind(dat_ecm, scrs_site) %>% 
  ggplot(aes(x=long, y=lat, colour=MDS1)) + 
  geom_point()


# including spatial patterns in analyses of community composition
# one way to do this is using principle coordinates of neighbour matrices
xy <- dist(dat_ecm[, c('long', 'lat')]) # create distance matrix of longs and lats
xy.pcnm <- pcnm(xy)$vectors %>% as.data.frame() # export the result into a dataframe

# combine the PCNM results with the rest of the data  
temp <- cbind(dat_ecm, xy.pcnm)

# what do the PCNMs represent? A plotting example
ggplot(temp, aes(x=long, y=lat, colour=PCNM1)) + 
  geom_point()

# what variation is explained by these spatial variables
cap.sp <- capscale(mat_ecm ~ ., data=temp %>% 
                     select('long', 'lat', starts_with('PCNM')))
cap.sp
# proportion of variation associated with each axis
cap.sp$CCA$eig/cap.sp$tot.chi
# explanatory value
anova(cap.sp)

# which individual spatial variables to include?
cap.0 <- capscale(mat_ecm ~ 1, data=xy.pcnm) # intercept-only, starting analysis
# uncomment this next line to run - takes a long time
# ordistep(cap.0, formula(cap.sp), direction='forward')


# should we include all climate variables in our analysis
# check for variance inflation - values higher than ~ 10 are unlikely to explain unique variation
cap.cl <- capscale(mat_ecm ~ tmean + prec + srad, data=dat_ecm)
vif.cca(cap.cl)

# variation partitioning - here to three groups of variables
vp <- varpart(vegdist(mat_ecm, distance='robust.aitchison'), 
              ~site_code, # type of environment sample was collected from = X1
              ~tmean + prec + srad, # climate = X2
              ~ long + lat + PCNM1 + PCNM2 + PCNM3 + PCNM10 + PCNM12 + # spatial = X3
                PCNM9 + PCNM13 + PCNM4, data=temp)
vp
plot(vp, Xnames=c('Site code', 'Climate', 'Space'))

# test unique variation explained by partitions using partial CAPs
# 1 - associated with site_code
anova(capscale(mat_ecm ~ site_code +
                 Condition(tmean + srad + prec + 
                             long + lat+ PCNM1 + PCNM2 + PCNM3 + PCNM10 + PCNM12 + PCNM9 + PCNM13 + PCNM4), data=temp, 
               distance='robust.aitchison'))
# 2 - associated with climate
anova(capscale(mat_ecm ~ tmean + srad + prec +
                 Condition(site_code + 
                             long + lat + PCNM1 + PCNM2 + PCNM3 + PCNM10 + PCNM12 + PCNM9 + PCNM13 + PCNM4), data=temp, 
               distance='robust.aitchison'))
# 3 - associated with space
anova(capscale(mat_ecm ~ long + lat + PCNM1 + PCNM2 + PCNM3 + PCNM10 + PCNM12 + PCNM9 + PCNM13 + PCNM4 +
                 Condition(tmean + srad + prec + site_code), data=temp, 
               distance='robust.aitchison'))




