
CNP_final<-read.csv("Stoich_CNP/CNP_final.csv")

library(ggplot2)
library(tidyverse)


CNP_long <- CNP_final %>%
  select(Sample, Fire.Interval, Fire.Severity, C_N, C_P, N_P) %>%
  filter(!Sample=='Summary') %>% 
  pivot_longer(cols = c(C_N, C_P, N_P), names_to = "Ratio", values_to = "Value") %>%
  mutate(
    Ratio = dplyr::recode(Ratio,
                   C_N = "C:N",
                   C_P = "C:P",
                   N_P = "N:P")
  )
library(ggplot2)
library(viridis)  # for visually distinct palettes (optional)

# Fire Interval Plot
Interval <- ggplot(CNP_long, aes(x = Fire.Interval, y = Value, color = Fire.Interval)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, linewidth=3) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 4) +
  facet_wrap(~ Ratio, scales = "free_y") +
  scale_color_manual(values = c("Short" = "#1f78b4", "Long" = "#33a02c")) +  # blue-green
  labs(
    colour = "Fire Frequency",  # Rename legend for color (Fire.Interval)
    x = "Fire Interval",
    y = "Ratio"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 36, color = 'black'),
    axis.text.y = element_text(size = 36, color = 'black'),
    axis.title.y = element_text(size = 36, color = 'black'),
    axis.title.x = element_text(size = 36),
    strip.text = element_text(size = 36, color = 'black'),  # <-- Add this
    legend.position = "right",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30)
  )

# Fire Severity Plot
Severity <- ggplot(CNP_long, aes(x = Fire.Severity, y = Value, color = Fire.Severity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, linewidth=3) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 4) +
  facet_wrap(~ Ratio, scales = "free_y") +
  scale_color_manual(values = c("Low" = "#e31a1c", "High" = "#ff7f00")) +  # red-orange
  labs(
    colour = "Fire Severity",  # Rename legend for color (Fire.Interval)
    x = "Fire Severity",
    y = "Ratio"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 36, color = 'black'),
    axis.text.y = element_text(size = 36, color = 'black'),
    axis.title.y = element_text(size = 36, color = 'black'),
    axis.title.x = element_text(size = 36),
    strip.text = element_text(size = 36, color = 'black'),  # <-- Add this
    legend.position = "right",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30)
  )

# Display both
Interval
Severity

library(patchwork)
p1<-Interval/Severity
p1

ggsave(filename = "plots/CNP_Regime.png", plot = p1, dpi=300, device = "png", width = 72, height = 40, units = "cm")
