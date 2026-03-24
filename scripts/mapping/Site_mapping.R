# ============================================================
# Publication-ready maps of experimental site locations
# ============================================================

# ---- Load required packages ----
library(tidyverse)
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# ---- 1. Read and clean site data ----

# (a) In-growth bag sites
Bag_sites <- read.csv("Processed_data/All_Bag_data.csv") %>%
  distinct(Site, Longitude, Latitude)

# (b) Second-round bag experiment
Bag_Site_2nd_rnd <- read_csv("ABS_Second_Rnd/Processed_data/Bag_Seq_wide.csv") %>%
  distinct(Site) %>%
  left_join(Bag_sites, by = "Site")

# (c) ABS main project sites
ABS_Site_Locations <- read_excel("Raw_data/Site_Data/60_Site_Locaitons.xlsx") %>%
  mutate(
    Site = str_remove(Site, "ABS0"),
    Site = as.numeric(Site)
  ) %>%
  distinct(Site, Longitude, Latitude)

# (d) Decomposition experiment sites
Decomp_Site_Locations <- read_csv("Raw_data/Site_Data/Decompisition Site Locations.csv") %>%
  rename(Site = `New Name`) %>%
  mutate(
    Site = str_remove(Site, "ABS0"),
    Site = as.numeric(Site)
  ) %>%
  distinct(Site, Longitude, Latitude)

# ---- 2. Get Australia base map ----
australia <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin == "Australia")

# ---- 3. Define plotting function ----
plot_sites <- function(data, title, shape, fill) {
  ggplot() +
    geom_sf(data = australia, fill = "grey95", color = "black", linewidth = 0.5) +
    geom_point(
      data = data,
      aes(x = Longitude, y = Latitude),
      shape = shape, size = 1.5, fill = fill, stroke = 1
    ) +
    coord_sf(
      xlim = c(150, 152), ylim = c(-32, -35.5),
      expand = FALSE
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey85", linetype = "dotted"),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title = element_blank(),
      axis.text = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(
      location = "tl", which_north = "true",
      height = unit(0.8, "cm"), width = unit(0.8, "cm")
    ) +
    ggtitle(title)
}

# ---- 4. Generate four maps ----
map_bag1 <- plot_sites(Bag_sites, "12 Sites", shape = 21, fill = "#1b9e77")
map_bag2 <- plot_sites(Bag_Site_2nd_rnd, "6 Sites", shape = 22, fill = "#d95f02")
map_abs <- plot_sites(ABS_Site_Locations, "63 Sites", shape = 23, fill = "#7570b3")
map_decomp <- plot_sites(Decomp_Site_Locations, "30 Sites", shape = 24, fill = "#e7298a")

map_bag1
map_bag2
map_abs
map_decomp
# ---- 5. Export each as high-resolution PNG ----
output_dir <- "Figures/Site_Maps"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "Map_Bag_Round1.png"), map_bag1, width = 6, height = 5, dpi = 600)
ggsave(file.path(output_dir, "Map_Bag_Round2.png"), map_bag2, width = 6, height = 5, dpi = 600)
ggsave(file.path(output_dir, "Map_ABS_Sites.png"), map_abs, width = 6, height = 5, dpi = 600)
ggsave(file.path(output_dir, "Map_Decomposition_Sites.png"), map_decomp, width = 6, height = 5, dpi = 600)

# ============================================================
# Each PNG is now ready to import into PowerPoint
# ============================================================
