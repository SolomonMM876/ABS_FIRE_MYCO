
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


#All meta data from 30 sites
Site_Precip_Temp_Elv <- read_excel("Processed_data/Site_Precip_Temp_Elv.xlsx")%>%
  rename(elev=wc2.1_30s_elev)%>%
  select(site,Annual_Temp,Annual_Prec,elev)%>%
  rename(Site=site)
Site_Precip_Temp_Elv$Site= sub(c('ABS00|ABS0'),'',Site_Precip_Temp_Elv$Site) 

Bag_Site<-read.csv('Processed_data/All_Bag_Site_Info.csv')


#Nute and Veg data
Nutrients_Transects<-read_excel('Processed_data/Nutrients_Transect_level.xlsx')
VEG_COVER_Transects <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                  sheet = "Transect.Level_Data")
VEG_COVER_Transects$Site= sub(c('ABS00|ABS0'),'',VEG_COVER_Transects$Site)
VEG_COVER_Transects$Transect= sub(c('T'),'',VEG_COVER_Transects$Transect)

Myco_plant_spp<-read.csv('Processed_data/Myco_host_abundance.csv')%>%
  mutate(Site=as.factor(Site),
         Transect=as.factor(Transect))


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
  left_join(VEG_COVER_Transects)%>%
  left_join(Nutrients_Transects)%>%
  select(-`Fire 3`,-`Interval (yrs)...16`,-`FESM severity category`)%>%
  left_join(Site_Precip_Temp_Elv)%>%
  left_join(Myco_plant_spp)%>%
  #this is selecting for the 30 decomp sites which I have CN data for all of
  filter(!is.na(Carbon))%>%
  #Site 60 was mislabeled and it wasnt clear what site 60 was so it is removed
  filter(!Site==60)


summary(Blast_ID)#no na;s in df



#funguild output
AM <- read_excel("Funguild/AM.xlsx")
Ecto <- read_excel("Funguild/Ecto.xlsx")
Path<- read_excel("Funguild/Path.xlsx")
Sap<- read_excel("Funguild/Sap.xlsx")
Myco<-read.csv("Funguild/Myco_guilds.csv")

#taxonomy table
otu<-read.csv('Raw_data/Updated_Data/demultiplexed.cleaned.combined.cf.fasta.blast.i97.a95.csv')
#remove first three rows that do not have taxonomy ids and 2 rows with totals and sample counts 
otu<-otu[-c(1:3),-c(2:3) ]

tax <- otu %>%
  separate(taxonomy, into = c("Percentage", "Species_Name", "Accession", "SH_ID", "Reference", "Taxonomy_Details"), sep = "\\|", extra = "drop") %>%
  separate(Taxonomy_Details, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  select(SH_ID:Species)%>%
  mutate(across(Kingdom:Species, ~ gsub('^[a-z]__', '', .)),
         guild = case_when(Genus %in% Myco$Genus ~ 'mycorrhizal',
                           #                Genus %in% Ecto$Genus ~ 'ectomycorrhizal', 
                           #Genus %in% AM$Genus ~ 'arbuscular_mycorrhizal', 
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
mat <- mat[rowSums(mat) >= 1000, ]  # Keep only samples with at least 1000 reads
temp <- rarecurve(mat, step=1000, tidy=TRUE)
Blast_ID %>%
  left_join(temp %>% rename(sample_ID = Site))%>%
  ggplot(aes(x=Sample, y=Species, colour=as.factor(Site), group=sample_ID)) + 
  geom_line()



# check variation in sample effort, looking for break points WHAT DOES THIS MEAN
hist(log10(rowSums(mat)))
sort(rowSums(mat))[1:63]





#alpha diversity measurements
myco_otus <- filter(tax, guild%in% c('mycorrhizal'))$SH_ID
mat_myco <-  filtered_mat %>% as.data.frame() %>% 
  dplyr::select(all_of(myco_otus))%>%# just ecto OTUs
  #remove samples that have no myco reads in them
  filter(!if_all(everything(), ~ . == 0))





alpha_diversity <- data.frame(
  Shannon = diversity(mat_myco, index = "shannon"),
  Simpson = diversity(mat_myco, index = "simpson"),
  Observed = rowSums(mat_myco > 0)  # Count of non-zero ASVs per sample
)

alpha_diversity <- mat_myco %>%
  rowwise() %>%  # Ensure calculations are row-wise
  mutate(
    Shannon = diversity(c_across(everything()), index = "shannon"),
    Simpson = diversity(c_across(everything()), index = "simpson"),
    Chao1 = chao1(c_across(everything())),
    Observed = sum(c_across(everything()) > 0),  # Count of non-zero ASVs
    Pielou = Shannon / log(Observed)  # Pielou's evenness
  ) %>%
  ungroup() %>%
  mutate(Pielou = ifelse(Observed <= 1, NA, Pielou)) %>%  # Avoid log(1) issue
  select(Shannon, Simpson, Chao1, Observed, Pielou)



alpha_diversity<- mat_myco%>%
  rownames_to_column('sample_ID')%>%
  select(sample_ID)%>%
  cbind(alpha_diversity)

alpha_myco_Site_Tran<-right_join(alpha_diversity,Blast_ID)%>%
  filter(!if_any(everything(), is.na))  # drop the sites that were removed because there was not enough seq depth



Shannon<-lmer((Shannon)~  Interval+Severity +Veg_Class + perc_myco_host_freq+
                NH4 + NO3 + Total.P +(1|Site) , #cant add Transect as a nested factor 
              data=alpha_myco_Site_Tran)

summary(Shannon)
Anova_resin<-round(Anova(Shannon,test='F'), 2) 
Anova_resin
plot(Shannon)
qqPlot(resid(Shannon))
#r2(Shannon)
emm_biomass_perc_myco_host_freq<-as.data.frame(emmeans(Shannon,
                                            ~perc_myco_host_freq))


library(glmmTMB)
Simpson <- lmer((Simpson) ~  Interval+Severity +Veg_Class +perc_myco_host_freq+
                     NH4 + NO3 + Total.P +(1|Site), 
                      data = alpha_myco_Site_Tran)


summary(Simpson)
Anova_resin<-round(Anova(Simpson,test='F'), 2) 
Anova_resin
plot(Simpson)
qqPlot(resid(Simpson))
r2(Simpson)
emm_biomass_Interval<-as.data.frame(emmeans(Simpson,
                                            ~Interval))

#Chao1
Chao1<-lmer(Chao1~  Interval+Severity +Veg_Class + perc_myco_host_freq+
              NH4 + NO3 + Total.P +(1|Site) , 
               data=alpha_myco_Site_Tran)


summary(Chao1)
Anova_resin<-round(Anova(Chao1,test='F'), 2) 
Anova_resin
plot(Chao1)
qqPlot(resid(Chao1))
r2(Chao1)
emm_biomass_Interval<-as.data.frame(emmeans(Chao1,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Chao1,
                                            ~Fire.Severity))
#Pielou
Pielou<-lmer(Pielou~  Interval+Severity +Veg_Class +perc_myco_host_freq+
              NH4 + NO3 + Total.P +(1|Site) , 
            data=alpha_myco_Site_Tran)


summary(Pielou)
Anova_resin<-round(Anova(Pielou,test='F'), 2) 
Anova_resin
plot(Pielou)
qqPlot(resid(Pielou))
r2(Pielou)
emm_biomass_Interval<-as.data.frame(emmeans(Pielou,
                                            ~perc_myco_host_freq))
emm_biomass_Severity<-as.data.frame(emmeans(Pielou,
                                            ~Fire.Severity))
