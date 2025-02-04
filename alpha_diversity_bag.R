
#LOAD LIBRARARIES
library(tidyr)
library(dplyr)
library(vegan)
library(tibble) 
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(fossil)


#All meta data from 12 sites with bags collected
Bag_data<-read.csv('Processed_data/All_Bag_Site_Info.csv')[,-1]%>%
  left_join(tmp<-read.csv('Processed_data/Resin_Nutrients.csv'))

Bag_data<-Bag_data%>%
  select(Tube_ID,Site, Transect, Location,Fire.Interval,Fire.Severity,#Site Data
         Ortho_P_mg_kg,Nitrate_mg_kg, Ammonia_mg_kg,#Nutrients from resin data
         Longitude,Latitude) #meta data

Seq_data<-read.csv('Processed_data/cleaned_seq_dat.csv')%>%
  #3 samples have no mycorrhizal fungi and are removed bringing total number to 88 samples
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))%>%
  filter(guild %in%c("Ectomycorrhizal","Arbuscular Mycorrhizal"))%>%
  select(-sample)

dat_summary_Jeff <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_summarised.tsv')%>%
  left_join(Seq_data%>%select(Tube_ID,barcode, Site, Transect, Location)%>%distinct())%>%filter(!is.na(Tube_ID))



Jeff_alpha<-Bag_data%>%left_join(dat_summary_Jeff, by= c('Tube_ID'))%>%select(Tube_ID,barcode)


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

Bag_Seq_wide<-left_join(Bag_data, #meta data
                        wide_seq)%>%
  mutate(across(starts_with('ITSall'), ~replace_na(., 0)))

#check to make sure there are no NAs
NA_rows <- Bag_Seq_wide %>%
  mutate(across(everything(), ~ ifelse(is.na(.), ., NA))) %>%
  filter(if_any(everything(), ~ !is.na(.)))

otu_table <-wide_seq%>%
  column_to_rownames(var='Tube_ID')#Remove SampleID column

all_zero_rows <- apply(otu_table, 1, function(x) all(x == 0))
which(all_zero_rows)  # Returns row indices where only one column has a nonzero value

single_nonzero_rows <- apply(otu_table, 1, function(x) sum(x > 0) == 1)
which(single_nonzero_rows)  # Returns row indices where only one column has a nonzero value


alpha_diversity <- otu_table %>%
  as.data.frame() %>%
  rowwise() %>%
  mutate(
    Shannon = diversity(c_across(everything()), index = "shannon"),
    Simpson = diversity(c_across(everything()), index = "simpson"),
    Chao1 = chao1(c_across(everything())),
    Observed = sum(c_across(everything()) > 0),  # Count of non-zero ASVs
    Pielou = Shannon / log(Observed)  # Pielou's evenness
  ) %>%
  ungroup() %>%
  mutate(Pielou = ifelse(Observed <= 1, NA, Pielou))%>%  # Avoid log(1) issue
select(Shannon,Simpson,Chao1,Observed,Pielou)

# Add SampleID back
alpha_diversity<-alpha_diversity%>%
  rownames_to_column(var='Tube_ID')%>%
  mutate(Tube_ID=as.numeric(Tube_ID))


# alpha_diversity<-alpha_diversity%>%
#   filter(!(Shannon  == 0 & Simpson  == 0))

tmp<-Bag_data%>%
  select(Tube_ID,Fire.Interval,Fire.Severity,Site, Transect, Location)

anti_join(tmp,wide_seq)

alpha_meta_Site_Tran_Loc<-left_join(Bag_data, alpha_diversity)

alpha_meta_Site<-left_join(Bag_data%>%select(-Transect,-Location), alpha_diversity)

alpha_meta_Site_Tran_Loc<-alpha_meta_Site_Tran_Loc%>%
  filter( !is.na(Shannon|Simpson))

#Shannon
Shannon<-lmer(scale(Shannon)~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                Nitrate_mg_kg+ Ammonia_mg_kg+ (1|Site/Transect) , 
                           data=alpha_meta_Site_Tran_Loc)


summary(Shannon)
Anova_resin<-round(Anova(Shannon,test='F'), 2) 
Anova_resin
plot(Shannon)
qqPlot(resid(Shannon))
#r2(Shannon)
emm_biomass_Interval<-as.data.frame(emmeans(Shannon,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Shannon,
                                            ~Fire.Severity))
#Pielous
Pielou<-lmer(sqrt(Pielou+0.01)~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                Nitrate_mg_kg+ Ammonia_mg_kg+ (1|Site/Transect) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Pielou)
Anova_resin<-round(Anova(Pielou,test='F'), 2) 
Anova_resin
plot(Pielou)
qqPlot(resid(Pielou))
#r2(Pielou)
emm_biomass_Interval<-as.data.frame(emmeans(Pielou,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Pielou,
                                            ~Fire.Severity))

Simpson<-lmer(sqrt(Simpson+0.01)~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                Nitrate_mg_kg+ Ammonia_mg_kg+ (1|Site/Transect) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Simpson)
Anova_resin<-round(Anova(Simpson,test='F'), 2) 
Anova_resin
plot(Simpson)
qqPlot(resid(Simpson))
r2(Simpson)
emm_biomass_Interval<-as.data.frame(emmeans(Simpson,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Simpson,
                                            ~Fire.Severity))

Chao1<-lmer(Chao1 ~  Fire.Severity+ Fire.Interval+Ortho_P_mg_kg+
                 Nitrate_mg_kg+ Ammonia_mg_kg+ (1|Site/Transect) , 
              data=alpha_meta_Site_Tran_Loc)


summary(Observed)
Anova_resin<-round(Anova(Observed,test='F'), 2) 
Anova_resin
plot(Observed)
qqPlot(resid(Observed))
r2(Observed)
emm_biomass_Interval<-as.data.frame(emmeans(Observed,
                                            ~Fire.Interval))
emm_biomass_Severity<-as.data.frame(emmeans(Observed,
                                            ~Fire.Severity))

