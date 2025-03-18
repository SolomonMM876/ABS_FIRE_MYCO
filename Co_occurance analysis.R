library(tidyverse)
library(cocorresp)
library(vegan)

source('Soil_ITS_Scripts/Soil_ITS_data prep.R')
#myco table
myco_com<-dat_myco%>%select(Site,Transect,SH_ID,readcount)%>%  unite(location, Site, Transect, sep = "_") %>% # Merge Site and Transect
pivot_wider(names_from = SH_ID, values_from = readcount, values_fill = 0)%>%column_to_rownames('location')
  
dat_myco%>%select(Site,Transect)%>%distinct()%>%unite(location, Site, Transect, sep = "_") %>%pull()->samps

meta<-Blast_ID%>%unite(location, Site, Transect, sep = "_",remove = FALSE) %>%select(Fire.Interval,Fire.Severity,location,Site,Transect)

#Veg data
ABS_Floristics<- read.csv("Raw_data/Site_Data/ABS_Floristic plot data_ANALYSES_V2_20240819_CEG.cleaned.csv")[,-1]%>%
  mutate(SiteNo = gsub('ABS00|ABS0', "", SiteNo))%>%
  rename(Site=SiteNo)

library(APCalign)
#load for austraits
tax_resources <- load_taxonomic_resources()

spp_list<-ABS_Floristics%>%distinct(ScientificName)

spp_list<-align_taxa(spp_list$ScientificName, resources = tax_resources)

upd_spp_list<-update_taxonomy(spp_list, taxonomic_splits = "most_likely_species",resources = tax_resources)%>%
  rename(taxon_name=`suggested_name`)%>%
  select(original_name,taxon_name,genus:taxon_rank)%>%
  filter(!is.na(taxon_name))


plant_comm <- ABS_Floristics %>%
  mutate(
    Transect1 = rowSums(select(., Q1.:Q15.), na.rm = TRUE),
    Transect2 = rowSums(select(., Q16.:Q30.), na.rm = TRUE)
  ) %>%  select(Site,ScientificName,Transect1,Transect2)%>%
  pivot_longer(cols = starts_with("Transect"), 
               names_to = "Transect", 
               values_to = "Freq") %>%
  mutate(Transect = gsub("Transect", "", Transect))%>%  # Remove "Transect" prefix
  unite(location, Site, Transect, sep = "_") %>% # Merge Site and Transect
  filter(location %in% samps)%>%
  group_by(location,ScientificName)%>%
  summarise(Freq=sum(Freq))%>%
  pivot_wider(names_from = ScientificName, values_from = Freq, values_fill = 0)%>%
  select(location, where(~ any(. != 0))) %>% # Remove columns with all zeros
  column_to_rownames('location')


#test

# Perform symmetric co-correspondence analysis
symmetric_coca <- coca(plant_comm ~.,data= myco_com, method = "symmetric")

# Extract eigenvalues from the Co-CA model
lambda_scores <- summary(symmetric_coca)$lambda  # Replace coca_model with your actual model object


lambda_1 <- round(lambda_scores[1], 3)  # Round for cleaner display
lambda_2 <- round(lambda_scores[2], 3)

#proportions <- round(lambda_scores / sum(lambda_scores) * 100, 2)


# Get axis correlations
#cor_axis <- corAxis(symmetric_coca)

# Permutation test for model significance
#perm_test <- permutest(symmetric_coca)

# Plot results
plot(symmetric_coca, main = "Symmetric Co-CA: Plant & Myco Fungi Communities")

# Extract site scores and reshape them
site_scores <- scores(symmetric_coca, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column(var = "location") %>%
  # Pivot longer to combine COCA.1 and COCA.2 for both Y and X
  pivot_longer(cols = starts_with("sites"),
               names_to = c("Type", "Axis"),
               names_pattern = "sites\\.([XY])\\.COCA\\.(\\d+)",
               values_to = "Score") %>%
  # Create the "Type" column with labels "plant" for Y and "myco" for X
  mutate(Type = ifelse(Type == "Y", "Plant", "Myco")) %>%
  # Optionally, reorder columns for clarity
  select(location, Type, Axis, Score) %>%
  pivot_wider(names_from = c(Axis), 
                values_from = Score)%>%
  rename(COCA_1 = `1`,COCA_2 = `2`)

# Extract species scores for Y
species_scores_Plants <- scores(symmetric_coca, display = "species")$species$Y %>%
  as.data.frame() %>%
  rownames_to_column(var = "Species") %>%
  mutate(Type = "Plant")

# Extract species scores for X
species_scores_Myco <- scores(symmetric_coca, display = "species")$species$X %>%
  as.data.frame() %>%
  rownames_to_column(var = "Species") %>%
  mutate(  Type = "Myco")

# Combine all scores into a single dataframe
species_scores <- bind_rows(species_scores_Plants, species_scores_Myco) %>%
  rename(COCA_1 = `COCA 1`, COCA_2 = `COCA 2`,Name=Species) %>%
  left_join(tax_myco%>%rename(Name=SH_ID)) %>%
  mutate(Name = if_else(!is.na(Genus), Genus, Name))

out_Interval
species_scores_Interval <- bind_rows(species_scores_Plants, species_scores_Myco) %>%
  rename(COCA_1 = `COCA 1`, COCA_2 = `COCA 2`,Name=Species) %>%
  left_join(out_Interval%>%rename(Name=SH_ID)) %>%
  mutate(Name = if_else(!is.na(Genus), Genus, Name))%>%
  filter(!is.na(Genus))

interval_colors <- c("Long" = "darkred", "Short" = "orange")

# first plot - site scores along with centroids for each group
right_join(meta, site_scores) %>% 
  ggplot(aes(x=COCA_1, y=COCA_2)) + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_text_repel(data=species_scores%>%
                    filter(abs(COCA_1) > 3.5 | abs(COCA_2) > 3.5),
                  aes(x=COCA_1, y=COCA_2, label=Name),
                  colour='black',size=4)+ 
  geom_point(data=species_scores,
               inherit.aes = FALSE,
               aes(x=COCA_1, y=COCA_2, color=Type), size=1, stroke = 3)+
  scale_shape_manual(values = c(3,1))+
  labs(x = paste0("COCA1 (λ = ", lambda_1, ")"),
       y = paste0("COCA2 (λ = ", lambda_2, ")")) +  # Insert raw lambda values
  theme_bw() + 
  theme(legend.position='top')->p2



p2















coca_passenger <- coca(plant_comm ~ ., data=myco_com, method = "predictive")


summary(coca_passenger)
 

coca_driver <- coca(myco_com ~., data= plant_comm, method = "predictive")

summary(coca_driver)

# Get percentage of variation explained
# Extract variance explained for COCA1 and COCA2
lambda1 <- coca_passenger$varianceExp$Yblock["Comp 1"]
lambda2 <- coca_passenger$varianceExp$Yblock["Comp 2"]

cat("Variance explained by COCA1:", lambda1, "%\n")
cat("Variance explained by COCA2:", lambda2, "%\n")

total_Y <- coca_passenger$totalVar$Yblock
total_X <- coca_passenger$totalVar$Xblock

cat("Total variance in Plant Community:", total_Y, "\n")
cat("Total variance in Mycorrhizal Community:", total_X, "\n")


# Cross-validation for AM fungi predicting plants
crossval_passenger <- crossval(plant_comm , myco_com)
crossval_driver <- crossval( myco_com,plant_comm)

print(crossval_passenger) # Should show proportion of variation predicted
print(crossval_driver)    # Compare to see which model is stronger

# Permutation test for predictive power
perm_test_passenger <- permutest(coca_passenger, permutations = 999)
perm_test_driver <- permutest(coca_driver, permutations = 999)

print(perm_test_passenger)
print(perm_test_driver)