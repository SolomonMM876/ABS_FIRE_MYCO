
#LOAD LIBRARARIES
library(tidyverse)
library(vegan)
library(tibble) 
library(readxl)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(fossil)



#veg datcar#veg data from sites
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
  left_join(Myco_plant_spp)%>%
  select(-`Fire 3`,-`Interval (yrs)...16`,-`FESM severity category`)

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
mat<-otu%>%
  select(-last_col())%>%
  remove_rownames()%>%
  column_to_rownames("SH_ID")%>%
  t()

# assess variation in sampling effort, plotting sample effort curves
filtered_mat <- mat[rowSums(mat) >= 3500, ]  # Keep only samples with at least 5000 reads

filtered_temp <- rarecurve(filtered_mat, step=1000, tidy=TRUE)

temp <- rarecurve(mat, step=1000, tidy=TRUE)

filter<-Blast_ID %>%
  left_join(filtered_temp %>% rename(sample_ID = Site))

unfilter<-Blast_ID %>%
  left_join(temp %>% rename(sample_ID = Site))

excluded<-anti_join(unfilter,filter)


depth_ID_filter<-filter%>%
  mutate(Site_Color = ifelse(as.numeric(as.character(Site)) %in% 49:63, "red", "black"),
         Site_Color = case_when(
           Site == 49 & Transect == 1 ~ "black",  # Change Site 49, Transect 1 to 'black'
           TRUE ~ Site_Color  # Keep all other values unchanged
         ))%>%
  ggplot(aes(x = Sample, y = Species, colour = Site_Color, group = sample_ID)) + 
  scale_color_identity() +  # Uses the colors directly without a legend
  geom_line() #+
  #theme(legend.position = "none")
depth_ID_filter

#plotly::ggplotly(depth_ID_filter)
depth_ID_excluded<-excluded%>%
  mutate(Site_Color = ifelse(as.numeric(as.character(Site)) %in% 50:63, "red", "black")) %>%
  ggplot(aes(x = Sample, y = Species, colour = Site_Color, group = sample_ID)) + 
  scale_color_identity() +  # Uses the colors directly without a legend
  geom_line()+
  theme(legend.position = "none")
depth_ID_excluded
plotly::ggplotly(depth_ID_filter)

library(patchwork)
depth_ID_excluded + depth_ID_filter


filter%>%select(Site,Transect)%>%distinct()

# check variation in sample effort, looking for break points WHAT DOES THIS MEAN
hist(log10(rowSums(filtered_mat)))
sort(rowSums(filtered_mat))[1:63]


#alpha diversity measurements
myco_otus <- filter(tax, guild%in% c('mycorrhizal'))$SH_ID
mat_myco <-  filtered_mat %>% as.data.frame() %>% 
  dplyr::select(all_of(myco_otus))



alpha_diversity <- mat_myco%>%
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


alpha_diversity<- alpha_diversity%>%
  #somehow row names are lost above, so I am taking them from original matrix and binding them back in
  cbind(mat_myco%>% rownames_to_column('sample_ID')%>%select(sample_ID))

alpha_myco_Site_Tran<-left_join(alpha_diversity,Blast_ID) %>%
  mutate(Transect=as.factor(Transect))




Shannon<-lmer((Shannon)~  Interval+Severity +Veg_Class +
                Tree.Basal.Area_m2 + Herb.Cover_0.50cm_perc + perc_myco_host_freq+(1|Site)  , 
              data=alpha_myco_Site_Tran)

Shannon<-lmer((Shannon)~  Interval+Severity +Veg_Class + perc_myco_host_freq+(1|Site)  , 
              data=alpha_myco_Site_Tran)

summary(Shannon)
Anova_resin<-round(Anova(Shannon,test='F'), 2) 
Anova_resin
plot(Shannon)
qqPlot(resid(Shannon))
#r2(Shannon)
emm_biomass_Interval_Shan<-as.data.frame(emmeans(Shannon,
                                            ~Interval))
emm_biomass_Interval_Shan

emm_biomass_Severity_Shan<-as.data.frame(emmeans(Shannon,
                                            ~Severity))

emm_biomass_Severity_Shan

# Simpson <- lmer(log10(Simpson+0.001) ~ Fire.Severity + Fire.Interval + (1|Site), 
#                       data = alpha_myco_Site_Tran)
# 
# library(car)
# summary(Simpson)
# Anova_resin<-round(Anova(Simpson,test='F'), 2) 
# Anova_resin
# plot(Simpson)
# qqPlot(resid(Simpson))
# r2(Simpson)
# emm_biomass_Interval<-as.data.frame(emmeans(Simpson,
#                                             ~Fire.Interval))
# emm_biomass_Severity<-as.data.frame(emmeans(Simpson,
#                                             ~Fire.Severity))



#Chao1
Chao1<-lmer(Chao1~  Interval+Severity +Veg_Class +
              Tree.Basal.Area_m2 + Herb.Cover_0.50cm_perc + perc_myco_host_freq+(1|Site) , 
               data=alpha_myco_Site_Tran)
Chao1<-lmer(Chao1~  Interval+Severity +Veg_Class + perc_myco_host_freq+(1|Site) , 
            data=alpha_myco_Site_Tran)

summary(Chao1)
Anova_resin<-round(Anova(Chao1,test='F'), 2) 
Anova_resin
plot(Chao1)
qqPlot(resid(Chao1))
emm_biomass_Interval_Chao1<-as.data.frame(emmeans(Chao1,
                                            ~Interval))
emm_biomass_Interval_Chao1

emm_biomass_Severity_Chao1<-as.data.frame(emmeans(Chao1,
                                            ~Severity))
emm_biomass_Severity_Chao1
#Pielou
Pielou<-lmer(Pielou~  Interval+Severity +Veg_Class +
               Tree.Basal.Area_m2 + Herb.Cover_0.50cm_perc + perc_myco_host_freq+(1|Site), 
             data=alpha_myco_Site_Tran)
Pielou<-lmer(Pielou~  Interval+Severity +Veg_Class + perc_myco_host_freq+(1|Site), 
             data=alpha_myco_Site_Tran)

summary(Pielou)
Anova_resin<-round(Anova(Pielou,test='F'), 2) 
Anova_resin
plot(Pielou)
qqPlot(resid(Pielou))
r2(Pielou)

emm_biomass_Severity<-as.data.frame(emmeans(Pielou,
                                            ~Severity))

