library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(readxl)

#OTU Coutibble#OTU Counts
SSU_NANOPORE <-read_table("Raw_data/Sequence/SSU/SSU_NANOPORE.csv",col_types = cols( .default = col_number()))
colnames(SSU_NANOPORE) <- gsub('"', '', colnames(SSU_NANOPORE))
SSU_NANOPORE$Site <- gsub(".*so", "", SSU_NANOPORE$Site)
SSU_NANOPORE<-SSU_NANOPORE[-3,] #I just removed one replicate so I can run the rest of this stuff

#Site Meta
Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')
Hyphal_CNHP <- read_excel("Raw_data/Stoich/Stoich_Totals_Round_1.xlsx")


Site_Nutrients<-Bag_Site%>%
  dplyr::select(Site,Bray.P:Fire.Interval,Tree.Basal.Area_m2:Dead.Tree.Canopy.Cover_perc,Pair)%>%
  dplyr::distinct()%>%
  group_by(Site) %>%
  summarise(across(c(Bray.P: Nitrogen), mean, na.rm = TRUE))





Site_Info<-read_excel("Raw_data/Site_Info_12_Short.xlsx")[,-1]%>%
  select(-Transect,-Sample)%>%distinct()%>% select(-`Fire 3`,- `FESM severity category`,
                                                   - `Interval (yrs)...13`, -`Interval (yrs)...14`)%>%
  left_join(Site_Nutrients)
  

#ID of OTU's
tax <- read_csv("Raw_data/Sequence/SSU/parsed_otus.csv")
  


mat<-SSU_NANOPORE%>%
  remove_rownames()%>%
  column_to_rownames("Site")

rarecurve(mat, step=1000, tidy=TRUE)%>%
  ggplot(aes(x=Sample, y=Species, group=Site)) + 
  geom_line() + 
  facet_wrap(~as.numeric(Site))

dat_ssu<-left_join(Site_Info,mat%>%rownames_to_column('Site'))%>%
  filter(!Site==56)
  


data_long <- pivot_longer(dat_ssu, cols = starts_with("OTU_"), names_to = "OTU_ID", values_to = "Count")%>%
  left_join(tax)


ggplot(data_long, aes(x = Site, y = Count, fill = genus)) +
  geom_bar(stat = 'identity', position = 'stack') +
  facet_wrap(~ Interval,scales = 'free_x') +
  theme_classic() +
  labs(title = "Species Counts at Each Site by Fire Interval",
       x = "Site",
       y = "Species Count")


# theme_classic()# check variation in sample effort, looking for break points
hist(log10(rowSums(mat)))
sort(rowSums(mat))[1:25]


# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from one or more groups of site_code
library(indicspecies)
multipatt(dat_ssu %>% select(starts_with('OTU_')), # first argument is the community table, select only those columns
          dat_ssu$Interval) -> res
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

adonis2(mat_ssu ~ 1, data=dat_ssu, distance='robust.aitchison', add=TRUE)

# try a different distance index - result is quite good
cap1 <- capscale(mat_ssu ~ Interval + Total.P, data=dat_ssu, distance='robust.aitchison', add=TRUE)
plot(cap1)
cap1 # summary of inertia
proportions_CAP<-round(cap1$CCA$eig/cap1$tot.chi *100, 1) # proportion of variation associated with each axis
proportions_MDS<-round(cap1$CA$eig/cap1$tot.chi *100, 1) # proportion of variation associated with each axis
anova(cap1)

cap1$tot.chi

scrs <- scores(cap1, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site <- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')
scrs_biplot <- scrs %>% filter(score=='biplot')

interval_colors <- c("Long" = "darkred", "Short" = "orange")


# first plot - site scores along with centroids for each group
p2<-cbind(dat_ssu, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2)) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_text(aes( label = label), color= 'black', size=3)+
  geom_point(aes( colour= Interval ,shape = Severity), size=8)+ 
  scale_colour_manual(values = interval_colors) +     # Custom colors for Interval
  #geom_text(data = scrs_cent, aes(label = label), size = 2) + 
 labs( x=  paste0("CAP1 (", proportions_CAP[1], "%)"), y=  paste0("CAP2 (", proportions_CAP[2], "%)"))+
  geom_segment(data=scrs_cent%>%
               filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5),
               inherit.aes = FALSE,
               aes(x=0,y=0, xend=CAP1, yend=CAP2, group=label),
               arrow = arrow(type = "closed",length=unit(3,'mm')),
               color= 'black') +
  geom_text_repel(data=scrs_cent%>%
                    filter(abs(CAP1) > 0.6 | abs(CAP2) > 0.6)%>%
                    rename(OTU_ID=label)%>%
                    left_join(tax),
                  inherit.aes = FALSE,
                  aes(x=CAP1, y=CAP2, label=OTU_ID),
                  colour='black',size=7)+

  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0.5,size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25) )+
  theme(legend.position='top')

p2
