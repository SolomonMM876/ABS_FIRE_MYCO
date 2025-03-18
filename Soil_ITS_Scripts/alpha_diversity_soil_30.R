source('Soil_ITS_Scripts/Soil_ITS_data prep.R')
library(fossil)
library(vegan)
library(performance)
library(lme4)
library(car)
library(emmeans)


#select only mycorrhizal OTUs from sites we previously selected
mat_myco <-  filtered_mat %>% as.data.frame() %>% 
  select(all_of(myco_otus))%>%# just ecto OTUs
  filter(rownames(.) %in% Site_Sample_IDs) # Keep only the desired samples



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


#extract rownames
alpha_diversity<- mat_myco%>%
  rownames_to_column('sample_ID')%>%
  select(sample_ID)%>%
  cbind(alpha_diversity)

alpha_myco_Site_Tran<-left_join(alpha_diversity,Blast_ID)

Shannon<-lmer((Shannon)~  Interval+Severity +(1|Site) , data=alpha_myco_Site_Tran)

summary(Shannon)
Anova_resin<-round(Anova(Shannon,test='F'), 2) 
Anova_resin
qqPlot(resid(Shannon))
plot(Shannon)
r2(Shannon)
emm_Shan_Interval<-as.data.frame(emmeans(Shannon,
                                            ~Interval))
emm_Shan_Interval

#simpson
hist(alpha_myco_Site_Tran$Simpson)
Simpson <- lmer((Simpson) ~  Interval+Severity+ +(1|Site), 
                      data = alpha_myco_Site_Tran)


summary(Simpson)
Anova_resin<-round(Anova(Simpson,test='F'), 2) 
Anova_resin
qqPlot(resid(Simpson))
plot(Simpson)
r2(Simpson)
emm_Simp_Interval<-as.data.frame(emmeans(Simpson,
                                            ~Interval))
emm_Simp_Interval

#Chao1
Chao1<-lmer(Chao1~  Interval+Severity +(1|Site) , 
               data=alpha_myco_Site_Tran)


summary(Chao1)
Anova_resin<-round(Anova(Chao1,test='F'), 2) 
Anova_resin
plot(Chao1)
qqPlot(resid(Chao1))
r2(Chao1)
emm_Chao_Interval<-as.data.frame(emmeans(Chao1,
                                            ~Interval))
emm_Chao_Interval
#Pielou
Pielou<-lmer(Pielou~  Interval+Severity +(1|Site) , 
            data=alpha_myco_Site_Tran)


summary(Pielou)
Anova_resin<-round(Anova(Pielou,test='F'), 2) 
Anova_resin
plot(Pielou)
qqPlot(resid(Pielou))
r2(Pielou)
emm_Pielou_Interval<-as.data.frame(emmeans(Pielou,
                                            ~Interval))
emm_Pielou_Interval

