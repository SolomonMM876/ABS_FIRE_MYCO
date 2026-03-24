library(dplyr)
library(readxl)
library(flextable)
library(officer)

# Load data
Nutrients.Site <- read_excel('Processed_data/Nutrients_Site_level.xlsx') # Soil nutrients
Nutrient_Resins <- read.csv('Processed_data/Resin_Nutrients.csv')        # Resin nutrients

# Summarise resin nutrients with SE
Resin_Sites <- Nutrient_Resins %>%
  group_by(Site,Fire.Interval,Fire.Severity) %>%
  summarise(
    PO4_Resin = mean(Ortho_P_mg_kg, na.rm = TRUE),
    PO4_SE = sd(Ortho_P_mg_kg, na.rm = TRUE) / sqrt(n()),
    NH4_Resin = mean(Ammonia_mg_kg, na.rm = TRUE),
    NH4_SE = sd(Ammonia_mg_kg, na.rm = TRUE) / sqrt(n()),
    NO3_Resin = mean(Nitrate_mg_kg, na.rm = TRUE),
    NO3_SE = sd(Nitrate_mg_kg, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Site = as.factor(Site))

# Combine with soil nutrient data
Nutrient_table <- left_join(Resin_Sites, Nutrients.Site, by = "Site")

# Format mean ± SE for resin nutrients (rounded: 1 decimal for mean, 1 for SE)
Nutrient_table_fmt <- Nutrient_table %>%
  mutate(
    `PO₄` = sprintf("%.1f ± %.1f", PO4_Resin, PO4_SE),
    `NH₄t` = sprintf("%.1f ± %.1f", NH4_Resin, NH4_SE),
    `NO₃t` = sprintf("%.1f ± %.1f", NO3_Resin, NO3_SE),
    Bray.P = round(Bray.P, 2),
    NH4 = round(NH4, 2),
    NO3 = round(NO3, 2),
    Total.P = round(Total.P, 2),
    Carbon = round(Carbon, 2),
  ) %>%
  # Select and rename final columns
  select(
    Site,
    Interval=Fire.Interval,Severity=Fire.Severity,
    `PO₄`,
    `NH₄t`,
    `NO₃t`,
    `Bray P` = Bray.P,
    `NH₄` = NH4,
    `NO₃` = NO3,
    `Total P` = Total.P,
    `Carbon`
  ) %>%
  mutate(Site=as.numeric(Site)) %>% 
  arrange(Site)

ft <- flextable(Nutrient_table_fmt) %>%
  set_header_labels(
    Site = " ",      # Blank column header
    `PO₄` = "PO₄",
    `NH₄t` = "NH₄",
    `NO₃t` = "NO₃",
    `Bray P` = "Bray P",
    `NH₄` = "NH₄",
    `NO₃` = "NO₃",
    `Total P` = "Total P",
    `Carbon` = "Carbon"
  ) %>%
  add_header_row(
    values = c("Site","","", "", "Resins (mg/kg)", "","","","Soil (mg/kg)", "", "")
  ) %>%
  vline(j =c(3,6), border = fp_border(color = "black", width = 1)) %>%
  autofit()

ft
# Save to Word
doc <- read_docx() %>%
  body_add_par("Table XXX. Nutrient concentrations from soil and resin extractions", style = "heading 2") %>%
  body_add_flextable(ft)

print(doc, target = "Tables/Nutrient_Table_Output.docx")
#########################################################


Ortho_P<- read.csv( "ABS_Second_Rnd/processed_data/Ortho_P_2nd_Rnd.csv")
Ortho_P <- head(Ortho_P, -2)

Ortho_P_Sec <- Ortho_P %>%
  select(Site, Transect, Location, Ortho_P_mg_kg) %>%
  mutate(across(c(Site, Transect, Location), as.factor))

# Calculate mean ± SE for each Site
Ortho_P_summary <- Ortho_P_Sec %>%
  group_by(Site) %>%
  summarise(
    Mean = mean(Ortho_P_mg_kg, na.rm = TRUE),
    SE = sd(Ortho_P_mg_kg, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(`Ortho-P` = sprintf("%.1f ± %.1f", Mean, SE)) %>%
  select(Site, `Ortho-P`) %>%
  arrange(as.numeric(as.character(Site)))  # Ensure proper numeric ordering

# Create the flextable
ft_ortho <- flextable(Ortho_P_summary) %>%
  add_header_row(values = c("Site", "Ortho-P (mg/kg)"), colwidths = c(1, 1)) %>%
  set_header_labels(Site = " ", `Ortho-P` = " ") %>%  # Blank column names
  autofit()

# Preview or save
ft_ortho
