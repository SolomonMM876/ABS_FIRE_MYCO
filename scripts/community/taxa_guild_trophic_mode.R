

tax<-read.csv('Processed_data/unknown_taxa_bags.csv')


library(fungaltraits)

fungal_traits<-fungal_traits()

myco_species_HMSC<-species_list%>%
  left_join(fungal_traits, by = 'species')


Ectomycorrhiza_lineage

Ectomycorrhiza_exploration_type

Growth_form




tax_with_traits <- tax %>%
  inner_join(fungal_traits, by = "species")

# Filter rows with no match in the first join
tax_unmatched <- tax_with_traits %>%
  filter(is.na(genus)) # Replace with a column that should not be NA if matched

# Extract letters before '_' in both dataframes
tax_unmatched <- tax_unmatched %>%
  mutate(prefix = sub("_.*", "", species))

fungal_traits <- fungal_traits %>%
  mutate(prefix = sub("_.*", "", species))

# Perform the second join by 'prefix'
tax_with_traits_updated <- tax_unmatched %>%
  left_join(fungal_traits, by = "prefix") %>%
  select(-prefix) %>%
  bind_rows(filter(tax_with_traits, !is.na(genus)))



# Install the PlutoF R package
install.packages("plutof")

# Load the package
library(plutof)

# Authenticate with PlutoF
plutof_auth()  # Follow the prompts to log in with your credentials
