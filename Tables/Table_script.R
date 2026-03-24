library(tidyverse)
library(knitr)
library(stringr)
library(officer)  # For Word output
library(flextable)

Bag_Perm<-read.csv('Tables/Bag_data_permanova.csv')
Soil_perm<-read.csv('Tables/Soil_permanova.csv')
Bag_Perm_2nd<-read.csv('Tables/Bag_data_2nd_permanova.csv')

# Clean factor names
factor_rename <- c(
  # "Ortho_P_mg_kg" = "Orthophosphate",
  # "Nitrate_mg_kg" = "Nitrate",
  # "Ammonia_mg_kg" = "Ammonia",
  "Fire.Interval" = "Fire Frequency",
  "Fire.Severity" = "Fire Severity",
  "Fire.Severity:Fire.Interval" = 'Fire Severity x Fire Frequency'
  
  # "All.Tree.Canopy.Cover_perc" = 'Canopy Cover',
  # "Tree.Basal.Area_m2" = "Tree Bsal Area",
  # "perc_myco_host_freq" = "Frequency of mycorrhizal host"
  
)

# Clean up and combine
perm_table <- bind_rows(Bag_Perm,Bag_Perm_2nd,Soil_perm)

perm_table_fmt <- perm_table %>%
  mutate(
    Factor = dplyr::recode(Factor, !!!factor_rename,.default = as.character(Factor)),
    p = ifelse(is.na(Pr..F.), NA, round(Pr..F., 3)),
    SumOfSqs = ifelse(is.na(SumOfSqs), NA, round(SumOfSqs, 1)),
    R2 = ifelse(is.na(R2), NA, round(R2, 2)),
    F = ifelse(is.na(F), NA, round(F, 2)),
  ) %>%
  select(Sample_Type, Factor, Df, SumOfSqs, R2, F, p)

# Build flextable
ft2 <- perm_table_fmt %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. perMANOVA results for fungal communities. Significant terms in bold (p ≤ 0.05).") %>%
  bold(i = ~ p <= 0.05, j = c("F", "p"), bold = TRUE) %>%
 # compose(i = ~ !duplicated(Sample_Type), j = "Sample_Type", value = as_paragraph(Sample_Type)) %>%
  hline(i = c(3,6), border = fp_border(width = 1), part = "body") %>%
  align(align = "left", part = "all") %>%
  autofit()

ft2

# Create a new Word document
doc <- read_docx()

# Add a title (optional)
doc <- doc %>%
  body_add_flextable(ft2)  # Replace `ft2` with your flextable object


print(doc, target = "Tables/perMANOVA_Table_Output.docx")


####Biomass table

Biomass_Aov<-read.csv('Tables/Biomass_Anova.csv')
Biomass__Second_Aov<-read.csv('Tables/Biomass_Anova_2nd_rnd.csv')

Biomass<-rbind(Biomass_Aov,Biomass__Second_Aov)

Biomass_table_fmt <- Biomass%>%
  mutate(
    Factor = dplyr::recode(Factor, !!!factor_rename,.default = as.character(Factor)),
    p = ifelse(is.na(Pr..F.), NA, round(Pr..F., 3)),
    F = ifelse(is.na(F), NA, round(F, 2)),
  ) %>%
  select(Round,Factor,F, Df, Df.res, p)

# Build flextable
ft3 <- Biomass_table_fmt %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Anova results for Biomass communities. Significant terms in bold (p ≤ 0.05).") %>%
  bold(i = ~ p <= 0.05, j = c("F", "p"), bold = TRUE) %>%
  compose(i = ~ !duplicated(Round), j = "Round", value = as_paragraph(Round)) %>%
  hline(i = 6, border = fp_border(width = 1), part = "body") %>%
  align(align = "left", part = "all") %>%
  autofit()

ft3

# Create a new Word document
doc <- read_docx()

# Add a title (optional)
doc <- doc %>%
  body_add_flextable(ft3)  # Replace `ft2` with your flextable object


print(doc, target = "Tables/Biomass_ANOVA_Table_Output.docx")



####Nutrient table

Nutri_Aov<-readRDS('Tables/Nutients_1st_Rnd_Anova.rds')
Nutri_Second_Aov<-readRDS('Tables/Nutients_2nd_Rnd_Anova.rds')

Nutrients<-rbind(Nutri_Aov,Nutri_Second_Aov)

Nutrient_table_fmt <- Nutrients%>%
  mutate(
    Factor = dplyr::recode(Factor, !!!factor_rename,.default = as.character(Factor)),
    p = ifelse(is.na(`Pr(>F)`), NA, round(`Pr(>F)`, 3)),
    F = ifelse(is.na(F), NA, round(F, 2)),
    Df.res = ifelse(is.na(Df.res), NA, round(Df.res, 1)),
  ) %>%
  select(Round,Nutrient,Factor,F, Df, Df.res, p) %>% 
  arrange(Round,Factor)

# Build flextable
ft4 <- Nutrient_table_fmt %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Anova results for available nutrients. Significant terms in bold (p ≤ 0.05).") %>%
  bold(i = ~ p <= 0.05, j = c("F", "p"), bold = TRUE) %>%
  #compose(i = ~ !duplicated(Round), j = "Round", value = as_paragraph(Round)) %>%
  hline(i = 6, border = fp_border(width = 1), part = "body") %>%
  align(align = "left", part = "all") %>%
  autofit()

ft4

# Create a new Word document
doc <- read_docx()

# Add a title (optional)
doc_nute <- doc %>%
  body_add_flextable(ft4)  # Replace `ft2` with your flextable object


print(doc_nute, target = "Tables/Nutrient_ANOVA_Table_Output.docx")



#Alpha diversity table#############
Anova_resin<-read.csv('Tables/Alpha_diversity_Anova.csv')
Anova_resin_2nd<-read.csv('Tables/Alpha_diversity_Anova_2nd_collection.csv')
Anova_soil<-read.csv('Tables/Alpha_diversity_Anova_soil.csv')


Anova_alpha<-bind_rows(Anova_resin,Anova_resin_2nd,Anova_soil)
Anova_alpha
Anova_resin_table_fmt <- Anova_alpha%>%
  mutate(
    Factor = dplyr::recode(Factor, !!!factor_rename,.default = as.character(Factor)),
    p = ifelse(is.na(`Pr..F.`), NA, round(`Pr..F.`, 3)),
    F = ifelse(is.na(F), NA, round(F, 2)),
    Df.res = ifelse(is.na(Df.res), NA, round(Df.res, 1)),
  ) %>%
  select(Source=source,Metric,Factor,F, Df, Df.res, p)

# Build flextable
ft5 <- Anova_resin_table_fmt %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Anova results for diversity metrics. Significant terms in bold (p ≤ 0.05).") %>%
  bold(i = ~ p <= 0.05, j = c("F", "p"), bold = TRUE) %>%
  compose(i = ~ !duplicated(Metric), j = "Metric", value = as_paragraph(Metric)) %>%
  hline(i = c(9,12), border = fp_border(width = 1), part = "body") %>%
  align(align = "left", part = "all") %>%
  autofit()

ft5

# Create a new Word document
doc <- read_docx()

# Add a title (optional)
doc_alpha_div <- doc %>%
  body_add_flextable(ft5)  # Replace `ft2` with your flextable object


print(doc_alpha_div, target = "Tables/Alpha_div_ANOVA_Table.docx")










