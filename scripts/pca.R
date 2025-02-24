library(tidyverse)
library(cowplot)
library(ggfortify)
library(ggrepel)
library(viridis)
library(ggbiplot)

setwd("C:/Users/13215/Box/Ch. 2 Intraspecific defense t-os")

locs = read_csv("data/localities.csv")

# pull climate data
climate = read_csv("data/flintbcm_climate_tall_herbarium.csv") %>% 
  filter(clim_year > 1950, clim_year < 2000) %>% 
  mutate(pck = abs(pck)) %>% 
  group_by(id) %>% 
  dplyr::summarize(cwd = sum(cwd), ppt_mm = sum(ppt_mm), pck = sum(pck), snw = sum(snw), str = sum(str), tmin = mean(tmin), tmax = mean(tmax)) %>% 
  left_join(., locs) %>% 
  filter(!is.na(cwd), taxon_name %in% c("Streptanthus tortuosus", "Streptanthus tortuosus var. tortuosus"))

table(climate$taxon_name)

# Define site groups
high_elev_sites <- c("LV3", "SQ3", "WL2", "YO1", "SQ1", "WL3", "SQ2", "YOSE8", "CP2", "YOSE10", "LV1", "LV2")
low_elev_sites <- c("WL1", "BH", "TM2", "KC2", "IH", "SC", "RB", "MtSH", "CALO", "TFC", "DRP", "LC")

# Function to calculate growing season length for high-elevation sites
calc_growing_season_high <- function(df) {
  df %>%
    arrange(clim_year, month) %>%
    group_by(id, clim_year) %>%
    summarize(
      start_month = first(month[pck == 0], default = NA),  # First month with no snowpack
      end_month = first(month[pck > 0 & !is.na(pck)], default = NA),  # First month with snow cover
      growing_season = ifelse(!is.na(start_month) & !is.na(end_month), end_month - start_month, NA)
    ) %>%
    ungroup()
}

# Function to calculate growing season length for low-elevation sites
calc_growing_season_low <- function(df) {
  df %>%
    arrange(clim_year, month) %>%
    group_by(id, clim_year) %>%
    summarize(
      start_month = first(month[ppt_mm > 25], default = NA),  # First month with ppt > 25 mm
      end_month = last(month[ppt_mm > 0], default = NA),  # Last month with any precipitation
      growing_season = ifelse(!is.na(start_month) & !is.na(end_month), end_month - start_month, NA)
    ) %>%
    ungroup()
}

climate_for_pc = climate %>% 
  select(cwd, pck, ppt_mm, snw, str, tmin, tmax, elevation) %>%

# Function to calculate growing season length for high-elevation sites
calc_growing_season_high <- function(df) {
  df %>%
    arrange(clim_year, month) %>%
    group_by(id, clim_year) %>%
    summarize(
      start_month = first(month[pck == 0], default = NA),  # First month with no snowpack
      end_month = first(month[pck > 0 & !is.na(pck)], default = NA),  # First month with snow cover
      growing_season = ifelse(!is.na(start_month) & !is.na(end_month), end_month - start_month, NA)
    ) %>%
    ungroup()
}
  
  
summary(climate_for_pc)

pc = prcomp(climate_for_pc, scale = TRUE)

pc_data = data.frame(pc$x)

locs_pc = cbind(climate, pc_data)

loadings = data.frame(varnames=rownames(pc$rotation), pc$rotation)

autoplot(pc, loadings = TRUE, loadings.label = TRUE, loadings.colour = "grey", loadings.label.colour = "black")
ggsave("figs/pca.pdf", width = 10, height = 10)

ggplot() +
  geom_point(data = filter(locs_pc, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = PC1, y = PC2), color = "purple") 
ggsave("figs/pca_with_sites.pdf", width = 10, height = 10)


ggplot() +
  geom_point(data = filter(locs_pc, type == "herbarium"), aes(x = PC1, y = PC2), color = "grey30", alpha = 0.3, size = 1) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = PC1, y = PC2), color = "purple") +
  geom_text_repel(data = filter(locs_pc, type == "seed_collection"), aes(x = PC1, y = PC2, label = id))
ggsave("figs/pca_with_sites_labeled.pdf", width = 10, height = 10)


ggplot() +
  geom_point(data = locs_pc, aes(x = PC1, y = PC2, color = -elevation), size = 1) +
  scale_color_viridis() +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = PC1, y = PC2), shape = 6) 
ggsave("figs/pca_with_elevation.pdf", width = 7, height = 7)

ggplot() +
  geom_point(data = locs_pc, aes(x = PC1, y = PC2, color = -latitude), size = 1, alpha = 0.5) +
  scale_color_viridis(option = "plasma") +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = PC1, y = PC2), shape = 6) 
ggsave("figs/pca_with_latitude.pdf", width = 7, height = 7)

ggplot() +
  geom_point(data = locs_pc, aes(y = tmin, x = elevation, color = -latitude), size = 1, alpha = 0.5) +
  scale_color_viridis(option = "plasma") +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(y = tmin, x = elevation), shape = 6) 
ggsave("figs/ppt_vs_tmin.pdf", width = 7, height = 7)

ggplot() +
  geom_point(data = locs_pc, aes(x = PC1, y = PC2, color = county)) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = PC1, y = PC2), shape = 6) 
ggsave("figs/pc_with_county.pdf", width = 7, height = 7)

# Low precipitation sites

low = filter(locs_pc, PC2 < -1)
write.csv(low, "data/low_ppt_locs.csv", row.names = FALSE)

high = filter(locs_pc, PC2 > 1)
write.csv(high, "data/high_ppt_locs.csv", row.names = FALSE)

high_el = filter(locs_pc, elevation > 3000)
write.csv(high_el, "data/high_el_locs.csv", row.names = FALSE)

snowy_site = filter(locs_pc, PC1 > 2)
write.csv(snowy_site, "data/snowy_site_locs.csv", row.names = FALSE)


# Make a map
states = map_data("state")

map = ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  coord_quickmap(xlim = c(-124, -118), ylim = c(36, 42)) +
  geom_point(data = filter(locs_pc, type == "herbarium"), aes(x = longitude, y = latitude), color = "grey", alpha = 0.5) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude), color = "purple", alpha = 1) +
  geom_text_repel(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude, label = id)) + 
  panel_border(size = 1, linetype = 1,
               remove = FALSE, colour = "black") +
  guides(colour = "none", fill = "none") +
  theme(panel.background = element_rect(fill = "grey", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.0, 0.1)); map
ggsave("figs/map.pdf", height = 10, width = 7)


map = ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  coord_quickmap(xlim = c(-124, -118), ylim = c(35, 42)) +
  geom_point(data = filter(locs_pc, type == "herbarium"), aes(x = longitude, y = latitude, color = PC1), size = 2) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude, color = PC1), size = 3) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude), shape = 1, size = 3) +
  # geom_text_repel(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude, label = id)) + 
  scale_color_viridis(option = "plasma", direction = -1) +
  panel_border(size = 1, linetype = 1,
               remove = FALSE, colour = "black") +
  theme(panel.background = element_rect(fill = "grey", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.01, 0.1)); map
ggsave("figs/map_PC1_nolabels.pdf", height = 10, width = 7)

map = ggplot() +
  geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "black", fill = "white") +
  coord_quickmap(xlim = c(-124, -118), ylim = c(35, 42)) +
  geom_point(data = filter(locs_pc, type == "herbarium"), aes(x = longitude, y = latitude, color = PC2), size = 2) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude, color = PC2), size = 3) +
  geom_point(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude), shape = 1, size = 3) +
  # geom_text_repel(data = filter(locs_pc, type == "seed_collection"), aes(x = longitude, y = latitude, label = id)) + 
  scale_color_viridis(option = "plasma", direction = -1) +
  panel_border(size = 1, linetype = 1,
               remove = FALSE, colour = "black") +
  theme(panel.background = element_rect(fill = "grey", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.01, 0.1)); map
ggsave("figs/map_PC2_nolabels.pdf", height = 10, width = 7)



elev = ggplot() +
  geom_point(data = filter(locs_pc, type == "herbarium", latitude > 35), aes(y = latitude, x = elevation), color = "grey", alpha = 0.5) + 
  geom_point(data = filter(locs_pc, type == "seed_collection", latitude > 35), aes(y = latitude, x = elevation), color = "purple") + 
  scale_fill_manual(values = c("turquoise", "coral")) +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = "none") +
  ylab("Latitude") +
  xlab("Elevation (m)"); elev
ggsave("figs/lat_vs_elevation", height = 6, width = 5)


