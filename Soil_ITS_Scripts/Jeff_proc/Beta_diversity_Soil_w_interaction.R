#LOAD LIBRARARIES
library(forcats)
library(tidyverse)
library(vegan)
library(indicspecies)
library(readxl)
library(ggplot2)



#From ITS prep script
Soil_Seq_wide<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_Seq_wide.csv')
myco_tax<-read.csv('Soil_ITS_scripts/Jeff_proc/Processed_data/Soil_dat_myco_tax.csv')

# next analysis - permanova

# first remove samples that have no mycorrhizal reads because they cause errors below
Soil_Seq_wide<-Soil_Seq_wide%>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

# extract the community table, save as a new object
mat_myco<- Soil_Seq_wide %>% select(starts_with("ITSall"))



temp<-adonis2(mat_myco~ Fire.Severity*Fire.Interval , data=Soil_Seq_wide, distance='robust.aitchison', by= 'margin')

temp<-temp%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>% 
  mutate(Sample_Type='Soil')
temp
write.csv(temp,'Tables/Soil_permanova.csv', row.names = FALSE)
table(Soil_Seq_wide$Fire.Interval)
table(Soil_Seq_wide$Fire.Severity)

#because there is a significant interaction I combining interval and severity into one factor and graphing this
Soil_Seq_wide <- Soil_Seq_wide %>%
  mutate(Fire = paste(Fire.Severity,"x", Fire.Interval, sep = "\n"))


cap.all <- capscale(mat_myco~ Fire
                    , data=Soil_Seq_wide, distance='robust.aitchison', add=TRUE)
anova(cap.all)
summary(cap.all)
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
cap.all
plot(cap.all)
proportions<-round(cap.all$CCA$eig/cap.all$tot.chi *100, 1) # proportion of variation associated with each axis

#cap_test <- capscale(mat_myco ~ 1 , data=Soil_Seq_wide, distance='robust.aitchison', add=TRUE)

#ordistep(cap_test, formula(cap.all), direction='forward')



# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap.all, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site<- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids') %>% 
  mutate(Fire=str_remove(label, 'Fire')) %>% 
  left_join(Soil_Seq_wide %>% distinct(Fire,Fire.Interval,Fire.Severity))%>%
  mutate(
    fill_color = case_when(
      Fire.Severity == "Low" ~ "white",
      Fire.Severity == "High" & Fire.Interval == "Long" ~ "darkred",
      Fire.Severity == "High" & Fire.Interval == "Short" ~ "orange"
    ),
    border_color = case_when(
      Fire.Interval == "Long" ~ "darkred",
      Fire.Interval == "Short" ~ "orange"
    )
  )
scrs_biplot <- scrs %>% filter(score=='biplot')





# regime_colors <- c(
#   "High\nx\nLong" = "darkred",
#   "High\nx\nShort" = "orange",
#   "Low\nx\nLong" = "white",
#   "Low\nx\nShort" = "white"
# )
# 
# legend_df <- tibble::tibble(
#   CAP1 = NA,
#   CAP2 = NA,
#   Fire = c("High\nx\nLong", "High\nx\nShort", "Low\nx\nLong", "Low\nx\nShort"),
#   fill_color = c("darkred", "orange", "white", "white"),
#   border_color = c("darkred", "orange", "darkred", "orange")
# )

# Ensure colors and linetypes are defined
interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_linetypes <- c("High" = "solid", "Low" = "dashed")

# Plot
cbind(Soil_Seq_wide, scrs_site) %>%
  ggplot(aes(x = CAP1, y = CAP2)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +  
  geom_point(aes(colour = Fire.Interval, shape = Fire.Severity), size = 10, stroke = 3) +
  stat_ellipse(aes(colour = Fire.Interval, linetype = Fire.Severity), size = 2) +
  labs(
    x = paste0("CAP1 (", proportions[1], "%)"), 
    y = paste0("CAP2 (", proportions[2], "%)"),
    colour = "Fire Frequency",
    linetype = "Fire Severity",
    shape = "Fire Severity"
  ) +
  scale_colour_manual(values = interval_colors) +
  scale_linetype_manual(values = severity_linetypes, guide='none') +
  scale_shape_manual(values = c("High" = 16, "Low" = 1)) +  # filled vs open points
  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 36, color = 'black'),
    axis.text.y = element_text(size = 36, color = 'black'),
    axis.title.y = element_text(size = 36, color = 'black'),
    axis.title.x = element_text(size = 36),
    legend.position = "top",
    legend.text = element_text(size = 30),        # Larger legend text
    legend.title = element_text(size = 30)        # Optional: larger legend title
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 8)),  # Adjust shape size
    colour = guide_legend(override.aes = list(size = 8))  # Adjust color legend size
  )  -> p2

p2

ggsave(filename = "plots/30_Soil_PCoA.png", plot = p2, dpi=300, device = "png", width = 65, height = 45, units = "cm")

