library(dplyr)
library(tidyr)


#taxonomy table
otu<-read.csv('Raw_data/Updated_Data/demultiplexed.cleaned.combined.cf.fasta.blast.i97.a95.csv')
#remove first three rows that do not have taxonomy ids and 2 rows with totals and sample counts 
otu<-otu[-c(1:3),-c(2:3) ]

tax <- otu %>%
  separate(taxonomy, into = c("Percentage", "Species_Name", "Accession", "SH_ID", "Reference", "Taxonomy_Details"), sep = "\\|", extra = "drop") %>%
  separate(Taxonomy_Details, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  select(SH_ID:Species)%>%
  mutate(across(Kingdom:Species, ~ gsub('^[a-z]__', '', .)))

# any kingdom-to-genus unclassified, change to NA
tax[tax %in% grep('unclassified', tax, value=T)] <- NA
# any species unclassified, change to NA
tax[tax$Species %in% grep('_sp', tax, value=T)] <- NA


library(stringr)

tax <- tax %>%
  mutate(Species = ifelse(str_detect(Species, '_sp'), NA, Species))


library(fungaltraits)
fungal_traits<-fungal_traits()

#na.omit(fungal_traits$em_expl)

duplicated_species <- tax$Species[duplicated(tax$Species)]



temp<-left_join(tax%>% rename(species=Species),
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