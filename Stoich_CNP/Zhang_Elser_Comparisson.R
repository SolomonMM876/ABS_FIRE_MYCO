library(lme4)
library(tidyverse)
library(car)
library(emmeans)
library(ggplot2)
library(readxl)

#samples 11.1 and 56.2 are weird, try analyses with and without

#run DF joins first
source('df_joins.R')

Bag_data<-read.csv('Processed_data/All_Bag_data.csv')

CNP_clean<-read.csv("Processed_data/CNP_clean.csv")

Zhang_Elser<-read_excel('Processed_data/Zhang_Elser_2017.xls')

#not currently incorporated
Panek_et_al_2024_CNP_mushrooms<-read_excel('Processed_data/Panek_et_al_2024_CNP_mushrooms.xlsx')


# Select relevant columns and pivot longer
CNP_clean_long <- CNP_clean %>%
  select(Sample, Carbon, Nitrogen, Percent_Phos_, C_N, C_P, N_P) %>%
  rename(P = Percent_Phos_,
         CN=C_N,
         C= Carbon,
         N=Nitrogen,
         CP=C_P,
         NP=N_P) %>%
  pivot_longer(cols = C:NP, names_to = "Element", values_to = "Value") %>%
  mutate(Source = "My Data")

Zhang_Elser_long_Ecto <- Zhang_Elser %>%filter(`Trophic Mode`=='Ectomycorrhizal')%>%
  select(C, N, P, CN, CP, NP) %>%
  pivot_longer(cols = C:NP, names_to = "Element", values_to = "Value") %>%
  mutate(Source = "Zhang Elser")

Zhang_Elser_long_All <- Zhang_Elser %>%
  select(C, N, P, CN, CP, NP) %>%
  pivot_longer(cols = C:NP, names_to = "Element", values_to = "Value") %>%
  mutate(Source = "Zhang Elser")

# Combine both datasets
CNP_comparison_Ecto <- bind_rows(CNP_clean_long, Zhang_Elser_long_Ecto)%>%
  filter(!is.na(Value))

CNP_comparison_All <- bind_rows(CNP_clean_long, Zhang_Elser_long_All)%>%
  filter(!is.na(Value))


#creat histograms
ggplot(CNP_comparison_All, aes(x = log10(Value))) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  facet_wrap(~Element, scales = "free") +  # Separate histogram for each element
  theme_minimal() +
  labs(title = "Distribution of Elemental Values",
       x = "Value",
       y = "Frequency")

comparison_results_All <- CNP_comparison_All %>%
  group_by(Element) %>%
  summarise(
    p_value =   t.test(log10(Value) ~ Source)$p.value) %>%
  mutate(label = paste0("p = ", signif(p_value, 3)))  # Show 3 significant digits

comparison_results_Ecto <- CNP_comparison_Ecto %>%
  group_by(Element) %>%
  summarise(
    p_value =   t.test(log10(Value) ~ Source)$p.value) %>%
  mutate(label = paste0("p = ", signif(p_value, 3)))  # Show 3 significant digits

# Create the boxplot
ggplot(CNP_comparison_Ecto, aes(x = Element, y = Value, fill = Source)) +
  geom_boxplot() +
  scale_y_log10() +  # Log scale for better visualization
  labs(title = "Comparison of Elemental Values with Ectos from Zhang and Elser 2017",
       x = "Element",
       y = "Value") +
  theme_minimal() +
  # Add p-values as text annotations
  geom_text(data = comparison_results_Ecto, aes(x = Element, y = max(CNP_comparison_Ecto$Value) * 1.1, label = label), 
            inherit.aes = FALSE, size = 5)

# Create the boxplot
ggplot(CNP_comparison_All, aes(x = Element, y = Value, fill = Source)) +
  geom_boxplot() +
  scale_y_log10() +  # Log scale for better visualization
  labs(title = "Comparison of Elemental with Ectos,Paths and Saprobs from Zhang and Elser 2017",
       x = "Element",
       y = "Value") +
  theme_minimal() +
  # Add p-values as text annotations
  geom_text(data = comparison_results_All, aes(x = Element, y = max(CNP_comparison_All$Value) * 1.1, label = label), 
            inherit.aes = FALSE, size = 5)
