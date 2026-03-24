library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(ggpubr)


#plot relationships between nuttrients

# Load necessary libraries
library(ggplot2)
library(ggpubr)

Bag_Site<-read.csv('Processed_data/All_Bag_data.csv')


# Define the variables of interest
vars <- c("Ammonia_mg_kg", "Nitrate_mg_kg", "Ortho_P_mg_kg", "pH")




# Fit model
m2 <- lmer(log10(Ammonia_mg_kg) ~ Fire.Severity + Fire.Interval + (1|Site/Transect), data = Bag_Site)

# Anova and summary
anova_res <- Anova(m2, test = "F")
summary(m2)

# Extract p-values for asterisks
p_interval <- anova_res["Fire.Interval", "Pr(>F)"]
p_severity <- anova_res["Fire.Severity", "Pr(>F)"]

# Get EMMs for Interval
emm_interval <- emmeans(m2, ~ Fire.Interval)
emm_df <- as.data.frame(emm_interval) %>%
  mutate(
    emmean_bt = 10^emmean,
    lower.CL_bt = 10^lower.CL,
    upper.CL_bt = 10^upper.CL
  )

# Add significance asterisk
sig_label <- case_when(
  p_interval < 0.001 ~ "***",
  p_interval < 0.01 ~ "**",
  p_interval < 0.05 ~ "*",
  TRUE ~ "ns"
)

# Plot
ggplot(emm_df, aes(x = Fire.Interval, y = emmean_bt)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower.CL_bt, ymax = upper.CL_bt), width = 0.2) +
  labs(
    title = "Nitrate vs Fire Interval",
    y = "Nitrate (mg/kg)",
    x = "Fire Interval"
  ) +
  annotate("text", x = 1.5, y = max(emm_df$upper.CL_bt) * 1.1, label = sig_label, size = 6) +
  theme_minimal()




# Create a correlation plot function
cor_plot <- function(x, y, data) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +  # Corrected column referencing
    geom_point(alpha = 0.6) + 
    geom_smooth(method = "lm", color = "blue", se = TRUE) + 
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + 
    theme_minimal() + 
    labs(title = paste("Correlation between", x, "and", y))
}

# Generate all pairwise correlation plots
plot_list <- list()
index <- 1
for (i in 1:(length(vars) - 1)) {
  for (j in (i + 1):length(vars)) {
    plot_list[[index]] <- cor_plot(vars[i], vars[j], Bag_data)
    index <- index + 1
  }
}

# Arrange plots together
ggarrange(plotlist = plot_list, ncol = 2, nrow = 3)





Nutrients_Transects<-read_excel('Processed_data/Nutrients_Transect_level.xlsx')
Bag_Site<-read.csv('Processed_data/All_Bag_data.csv')

Twelve_sites <- semi_join( Nutrients_Transects,Bag_Site%>%
                             mutate(Transect=as.factor(Transect),
                                    Site=as.factor(Site))
                           
                           , by = c("Transect", "Site"))%>%
  arrange(as.numeric(Site))%>%
  left_join(Bag_Site%>%select(Site,Transect,Fire.Interval)%>%
              mutate(Transect=as.factor(Transect),
                     Site=as.factor(Site)))

Twelve_sites$Site<-as.factor(Twelve_sites$Site)


 Twelve_long<- Twelve_sites %>%
  pivot_longer(
    cols = Bray.P:Nitrogen,
    names_to = "Nutrient",
    values_to = "Value" )%>%
  arrange(as.numeric(Site))
 
facet_order <- c("NO3", "NH4", "Nitrogen", "Bray.P", "Total.P", "Carbon")
 
 
Twelve_long%>%
ggplot( aes(x = Site, y = Value)) +
  geom_boxplot() +
  facet_wrap(~ Nutrient, scales = "free_y") +
  labs(x = "Site", y = "Nutrient Value (mg/kg)") +
  theme_minimal()

# Arrange the data by Site
Twelve_sites <- Twelve_sites %>%
  mutate(Site = factor(Site, levels = c("5", "7", "8", "10", "11", "12", "26", "29", "31", "34", "49", "56")))

# Create each plot individually
plot_NO3 <- ggplot(Twelve_sites, aes(x = Site, y = NO3)) +
  geom_boxplot() +
  labs(x = NULL, y = "NO3 (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_NH4 <- ggplot(Twelve_sites, aes(x = Site, y = NH4)) +
  geom_boxplot() +
  labs(x = NULL, y = "NH4 (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_Nitrogen <- ggplot(Twelve_sites, aes(x = Site, y = Nitrogen)) +
  geom_boxplot() +
  labs(x = NULL, y = "Nitrogen (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_BrayP <- ggplot(Twelve_sites, aes(x = Site, y = Bray.P)) +
  geom_boxplot() +
  labs(x = NULL, y = "Bray.P (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_TotalP <- ggplot(Twelve_sites, aes(x = Site, y = Total.P)) +
  geom_boxplot() +
  labs(x = "Site", y = "Total.P (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

plot_Carbon <- ggplot(Twelve_sites, aes(x = Site, y = Carbon)) +
  geom_boxplot() +
  labs(x = "Site", y = "Carbon (mg/kg)") +
  facet_wrap(~ Fire.Interval, scales = "free_x") +
  theme_minimal()

# Arrange plots using ggarrange and remove specific elements
ggarrange(
  plot_NO3 + rremove("x.text") + rremove("x.axis"),
  plot_NH4 + rremove("x.text") + rremove("x.axis"),
  plot_Nitrogen + rremove("x.text") + rremove("x.axis"),
  plot_BrayP + rremove("x.text") + rremove("x.axis"),
  plot_TotalP + rremove("y.text") + rremove("x.axis"),
  plot_Carbon + rremove("y.text") + rremove("y.axis"),
  ncol = 2, nrow = 3
)
