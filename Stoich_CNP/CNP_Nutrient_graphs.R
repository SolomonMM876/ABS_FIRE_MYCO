library(ggplot2)
library(lme4)
library(broom.mixed)
library(dplyr)
library(tidyr)
library(tibble)
library(MuMIn)
library(lmerTest)  # Enables F-tests for lmer models
library(car)

#run DF joins first
source('df_joins.R')

Bag_data<-read.csv('Processed_data/All_Bag_data.csv')

tmp<-Bag_data %>% 
  filter((Site==56 &Transect==2))

CNP_clean<-read.csv("Processed_data/CNP_clean.csv")

Bag_data<-Bag_data%>%
  group_by(Site, Transect) %>%
  summarise(Ortho_P_mg_kg = mean(Ortho_P_mg_kg, na.rm = TRUE),
            Nitrate_mg_kg = mean(Nitrate_mg_kg, na.rm = TRUE),
            Ammonia_mg_kg = mean(Ammonia_mg_kg, na.rm = TRUE)) %>% 
  left_join(CNP_clean)




# Define measured response variables
response_vars <- c("Carbon", "Nitrogen", "Percent_Phos_", "C_N", "C_P", "N_P")

# Gather data into long format
Bag_long <- Bag_data %>%
  pivot_longer(cols = all_of(response_vars), names_to = "Response", values_to = "Value")

# Define response variables and predictors
response_vars <- c("Carbon", "Nitrogen", "Percent_Phos_", "C_N", "C_P", "N_P")
predictors <- c("Ortho_P_mg_kg", "Nitrate_mg_kg", "Ammonia_mg_kg", 
                "Fire.Interval", "Fire.Severity")

# Initialize list to store results
result_list <- list()

# Loop over each response variable
for (resp in response_vars) {
  
  # Prepare data
  dat <- Bag_data %>%
    select(all_of(c("Site", "Ortho_P_mg_kg", "Nitrate_mg_kg", "Ammonia_mg_kg",
                    "Fire.Interval", "Fire.Severity", resp))) %>%
    rename(Value = all_of(resp)) %>%
    mutate(Value = log10(Value)) %>%
    filter(!is.na(Value))
  
  # Fit model
  model <- lmer(Value ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + 
                  Fire.Interval + Fire.Severity + (1|Site), data = dat)
  
  # Extract ANOVA table
  anova_res <- Anova(model, test='F') %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term %in% predictors) %>%
    mutate(Response = resp,
           F_value = round(`F`, 2),
           P_value = round(`Pr(>F)`, 3),
           Df_numerator = Df,
           Df_residual = Df.res )%>%
    select(Response, term, Df_numerator, Df_residual, F_value, P_value)
  
  result_list[[resp]] <- anova_res
}

# Combine all results into one tidy dataframe
tidy_results <- bind_rows(result_list)

# View table
print(tidy_results)

library(flextable)

# Clean factor names
factor_rename <- c(
  "Ortho_P_mg_kg" = "Orthophosphate",
  "Nitrate_mg_kg" = "Nitrate",
  "Ammonia_mg_kg" = "Ammonia",
  "Fire.Interval" = "Fire frequency",
  "Fire.Severity" = "Fire severity"
)

# Clean response names
response_rename <- c(
  "Carbon" = "% Carbon",
  "Nitrogen" = "% Nitrogen",
  "Percent_Phos_" = "% Phosphorous",
  "C_N" = "C:N",
  "C_P" = "C:P",
  "N_P" = "N:P"
)

linear_model_table <- tidy_results %>%
  mutate(
    Response = ifelse(Response %in% names(response_rename), response_rename[Response], Response),
    Factor = ifelse(term %in% names(factor_rename), factor_rename[term], term),
    p = round(P_value, 2),
    Df = round(Df_numerator, 2),
    Df.res = round(Df_residual, 1),
    F = round(F_value, 2)
  ) %>%
  select(Response, Factor, F,Df,Df.res, p)

# Build flextable
ft1 <- linear_model_table %>%
  flextable() %>%
  theme_booktabs() %>%
  bold(i = ~ p <= 0.05, j = c("F", "p"), bold = TRUE) %>%
  add_header_lines("Table X. Linear mixed model results. Significant terms in bold (p ≤ 0.05).") %>%
  compose(i = ~ !duplicated(Response), j = "Response", value = as_paragraph(Response)) %>%
  hline(i = which(!duplicated(linear_model_table$Response))+4, part = "body") %>%
  align(align = "left", part = "all") %>%
  autofit()

ft1


library(officer)
# Create a new Word document
doc <- read_docx()

# Add a title (optional)
doc <- doc %>%
  body_add_flextable(ft1)  # Replace `ft2` with your flextable object


print(doc, target = "Tables/CNP_Sig_Table_Output.docx")





