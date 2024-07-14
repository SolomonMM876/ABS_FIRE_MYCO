# install.packages("devtools")
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits)
fungal_traits<-fungal_traits()



na.omit(fungal_traits$em_expl)


temp<-left_join(out%>% rename(species=Species),
          fungal_traits%>%select(species,Genus,em_expl,em_text), by= 'species')
out$Genus
fungal_traits$species


# Perform a join based on the species_prefix
joined_by_species <- out %>%
  rename(species=Species) %>%
  left_join(fungal_traits %>% select(species, Genus, em_expl, em_text), by = 'species')

# Perform a join based on the Genus column
joined_by_genus <- out %>%
  rename(species=Species) %>%
  left_join(fungal_traits %>% select(species, Genus, em_expl, em_text), by = 'Genus')
 

tmp<-left_join( dat_ecm_12_site %>% 
  pivot_longer(cols = ends_with(".09FU"), names_to = "SH_ID", values_to = "value") 