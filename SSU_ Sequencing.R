library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(readxl)

#OTU Coutibble#OTU Counts
SSU_NANOPORE <-read_table("Raw_data/Sequence/SSU_NANOPORE.csv",col_types = cols( .default = col_number()))
colnames(SSU_NANOPORE) <- gsub('"', '', colnames(SSU_NANOPORE))
SSU_NANOPORE$Site <- gsub(".*so", "", SSU_NANOPORE$Site)
SSU_NANOPORE<-SSU_NANOPORE[-3,] #I just removed one replicate so I can run the rest of this stuff

#Site Meta
Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')
Hyphal_CNHP <- read_excel("Raw_data/Stoich/Stoich_Totals_Round_1.xlsx")


Bag_Site_Short<-Bag_Site%>%
  mutate(Regime = paste(Fire.Interval, Fire.Severity, sep = "_"))%>%
  dplyr::select(Site,Transect,Bray.P:Dead.Tree.Canopy.Cover_perc,Pair,log10_Second_Weight_bag_yield_est:Regime)%>%
  group_by(Site,Transect)%>%
  mutate(Biomass= mean(.085+10^log10_Second_Weight_bag_yield_est, na.rm = TRUE))%>%
  select(-log10_Second_Weight_bag_yield_est)%>%distinct()%>%
  left_join(Hyphal_CNHP%>%rename(Carb_Hyph = Carbon, Nitrog_Hyph =Nitrogen, Phos_Hyph= Percent_Phos)%>%
              select(-Fire.Interval,-Fire.Severity))





Site_Info<-read_excel("Raw_data/Site_Info_12_Short.xlsx")[,-1]%>%
  select(-Transect,-Sample)%>%distinct()%>% select(-`Fire 3`,- `FESM severity category`,
                                                   - `Interval (yrs)...13`, -`Interval (yrs)...14`)
  

#ID of OTU's
parsed_otus <- read_csv("Raw_data/Sequence/parsed_otus.csv")
  


mat<-SSU_NANOPORE%>%
  remove_rownames()%>%
  column_to_rownames("Site")

rarecurve(mat, step=1000, tidy=TRUE)%>%
  ggplot(aes(x=Sample, y=Species, group=Site)) + 
  geom_line() + 
  facet_wrap(~as.numeric(Site))

dat_ssu<-left_join(Site_Info,mat%>%rownames_to_column('Site'))%>%
  filter(!Site==56)
  

# check variation in sample effort, looking for break points
hist(log10(rowSums(mat)))
sort(rowSums(mat))[1:25]


# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from one or more groups of site_code
library(indicspecies)
multipatt(dat_ssu %>% select(starts_with('OTU_')), # first argument is the community table, select only those columns
          dat_ssu$Site) -> res
summary(res)

# to visualise differences in taxonomic composition, using output from multipatt()
# first prepare the data in the res object so that it can be joined with taxonomic info
# output is only those otus significant for one or more site_codes
out <- res[['sign']] %>% 
  filter(p.value <= 0.05) %>% 
  rownames_to_column('OTU_ID') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(site_code = gsub('^s.', '', group))

# next analysis - permanova
# extract the community table, save as a new object
mat_ssu <- dat_ssu %>% select(starts_with('OTU_'))

adonis2(mat_ssu ~ Site, data=dat_ssu, distance='robust.aitchison', add=TRUE)

# try a different distance index - result is quite good
cap1 <- capscale(mat_ssu ~ Site, data=dat_ssu, distance='robust.aitchison', add=TRUE)
plot(cap1)
cap1 # summary of inertia
cap1$CCA$eig/cap1$tot.chi # proportion of variation associated with each axis
anova(cap1)

