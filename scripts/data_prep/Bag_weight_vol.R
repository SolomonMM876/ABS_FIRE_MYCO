library(readxl)
library(tidyr)
library(dplyr)

library(ggplot2)


Fresh_Bag_weights <- read_excel("C:/Users/90957135/OneDrive - Western Sydney University/Bag_data/Fresh_Bag_weights.xlsx")

Dry_Bag_weights <- read_excel("C:/Users/90957135/OneDrive - Western Sydney University/Bag_data/empty_bag_w.xlsx",  2)

Dry_Bag_weights%>%
  summarise(
    mean_value = mean(wt_g, na.rm = TRUE),
    se_value = sd(wt_g, na.rm = TRUE) / sqrt(n())
  )

Fresh_Bag_weights<-Fresh_Bag_weights%>%
  group_by(Condition)%>%
  mutate(mean_weight= mean(weight_g, na.rm=T),
         std_err_weight= sd(weight_g, na.rm = TRUE) / sqrt(n()),
         mean_vol = mean(Volume,na.rm = T),
         std_err_vol= sd(Volume, na.rm = TRUE) / sqrt(n()),
         percentage_error_weight = (std_err_weight / mean_weight) * 100,
         percentage_error_vol = (std_err_vol / mean_vol) * 100)
#moisture_content
 Fresh_Bag_weights %>%
  group_by(Condition) %>%
  summarize(mean_weight = mean(mean_weight, na.rm = TRUE),) %>%
  pivot_wider(names_from = Condition, values_from = mean_weight) %>%
  mutate(moisture_content = ((Normal - Dry) / Dry) * 100) %>%
  pull(moisture_content)

moisture_content <- Fresh_Bag_weights %>%
  group_by(Condition) %>%
  summarize(mean_weight = mean(mean_weight, na.rm = TRUE),) %>%
  pivot_wider(names_from = Condition, values_from = mean_weight) %>%
  mutate(moisture_content = ((Normal - Dry) / Dry) * 100) %>%
  pull(moisture_content)


Fresh_Bag_weights%>%
  filter(Condition=='Dry')%>%
  summarise(dry_w= mean(weight_g),
            se_value = sd(weight_g, na.rm = TRUE) / sqrt(n())
  )



# Prepare Dry_Bag_weights
dry_bag_summary <- Dry_Bag_weights %>%
  summarise(
    weight = mean(wt_g, na.rm = TRUE),
    se = sd(wt_g, na.rm = TRUE) / sqrt(n())
  ) %>%
  mutate(source = "Dry Bag")

# Prepare Fresh_Bag_weights (filter for Condition == 'Dry')
fresh_bag_summary <- Fresh_Bag_weights %>%
  filter(Condition == "Dry") %>%
  summarise(
    weight = mean(weight_g, na.rm = TRUE),
    se = sd(weight_g, na.rm = TRUE) / sqrt(n())
  ) %>%
  mutate(source = "Fresh Bag")

# Combine the two data frames
combined_weights <- bind_rows(dry_bag_summary, fresh_bag_summary)

# Create the box plot
ggplot(combined_weights, aes(x = source, y = weight, fill = source)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = weight - se, ymax = weight + se), width = 0.2) +
  labs(
    title = "Comparison of Weights Between Dry and Fresh Bags",
    x = "Source",
    y = "Weight (g)"
  ) +
  theme_minimal()








# Box plots for the two bag conditions
 ggplot(Fresh_Bag_weights, aes(x = Condition, y = weight_g)) +
  geom_boxplot() +
  labs(x = "Condition", y = "Weight (g)", title = "Box Plot of Bag Weights")+
   theme_minimal()

# Bar graph with standard error
ggplot(Fresh_Bag_weights, aes(x = Condition, y=weight_g)) +
  geom_bar(stat="summary")+
  geom_errorbar(aes(ymin = mean_weight - std_err_weight, ymax = mean_weight + std_err_weight),
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Condition", y = "Mean Weight (g)", title = "Bar Graph of Bag Conditions with Standard Error")+
  theme_minimal()


Fresh_Bag_weights%>%
  ggplot(aes(x=Condition,y=weight_g))+
  geom_bar(stat="summary")+
  geom_errorbar(aes(ymin = mean_weight - std_err_weight, ymax = mean_weight + std_err_weight),
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Condition", y = "Mean Volumes cm^3", title = "Bar Graph of Bag Conditions with Standard Error")+
  theme_minimal()
