library(dplyr)
library(tidyverse)
library(flextable)
library(officer)

# Load data
CNP_clean <- read.csv('Stoich_CNP/CNP_final.csv')

# Ensure proper column names
CNP_clean <- CNP_clean %>%
  rename(
    C = Carbon,
    N = Nitrogen,
    P = Percent_Phos_,
    Frequency = Fire.Interval,
    Severity = Fire.Severity
  )

# Define helper function to calculate mean ± SE
mean_se <- function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  sprintf("%.1f ± %.1f", m, se)
}

mean_se_sci <- function(x) {
  m <- signif(mean(x, na.rm = TRUE), 1)
  se <- signif(sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))), 1)
  paste0(format(m, scientific = TRUE), " ± ", format(se, scientific = TRUE))
}
grouped_data <- bind_rows(
  CNP_clean %>% mutate(Grouping = "Overall"),
  CNP_clean %>% filter(Frequency == "Short") %>% mutate(Grouping = "Short Interval"),
  CNP_clean %>% filter(Frequency == "Long") %>% mutate(Grouping = "Long Interval"),
  CNP_clean %>% filter(Severity == "High") %>% mutate(Grouping = "High Severity"),
  CNP_clean %>% filter(Severity == "Low") %>% mutate(Grouping = "Low Severity")
)

#convert total P to ppm
grouped_data$P_ppm <- grouped_data$P * 10000


# Summarise for each group
summary_table <- grouped_data %>%
  group_by(Grouping) %>%
  summarise(
    `% C` = mean_se(C),
    `% N` = mean_se(N),
    `P ppm` = mean_se(P_ppm),
    `C:N` = mean_se(C_N),
    `C:P` = mean_se(C_P),
    `N:P` = mean_se(N_P),
    .groups = "drop"
  ) %>%
  mutate(Grouping = factor(Grouping, levels = c(
    "Overall", "Short Frequency", "Long Frequency", "High Severity", "Low Severity"
  ))) %>%
  arrange(Grouping)

# Convert to flextable
ft <- summary_table %>%
  flextable() %>%
  autofit()

ft
# Save to Word
doc <- read_docx() %>%
  body_add_par("Table XXX Elemental composition of mycorrhizal fungal biomass collected", style = "heading 2") %>%
  body_add_flextable(ft)

print(doc, target = "Tables/CNP_Table_Output.docx")
