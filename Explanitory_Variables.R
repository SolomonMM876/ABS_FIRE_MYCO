

library(readxl)
library(dplyr)
library(car)
library(lme4)

Bag_Site<-read_excel('Processed_data/All_Bag_Site_Info.xlsx')

Bag_Site_Short<-Bag_Site%>%
  mutate(Regime = paste(Fire.Interval, Fire.Severity, sep = "_"))%>%
  dplyr::select(Bray.P:Nitrogen,Tree.Basal.Area_m2:Dead.Tree.Canopy.Cover_perc,log10_Second_Weight_bag_yield_est)


explanitory_variables<-colnames(Bag_Site_Short)

# Perform ANOVA for each nutrient and store the results
anova_results <- lapply(explanitory_variables, function(nutrient) {
  formula <- as.formula(paste(nutrient, "~ Fire.Interval + Fire.Severity+ (1|Site)"))
  model <- lme4::lmer(formula, data = Bag_Site)
  Anova(model, test = "F")
})

# Extract p-values for Interval and Severity
p_values <- sapply(anova_results, function(result) {
  anova_table <- result
  c(Interval = anova_table$`Pr(>F)`[1],
    Severity = anova_table$`Pr(>F)`[2])
})


# Format p-values into a data frame
p_values_df <- as.data.frame(t(p_values))
colnames(p_values_df) <- c("Interval", "Severity")
rownames(p_values_df) <- explanitory_variables

# Print the resulting table
print(p_values_df)

p_values_df<-p_values_df%>%
  rownames_to_column(var='Variable')


library(writexl)

write_xlsx(p_values_df,'Processed_data/Explan_var_P_Regime.xlsx')
