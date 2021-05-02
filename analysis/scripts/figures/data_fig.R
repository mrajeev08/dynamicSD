# Data figure ------

library(data.table)
library(dplyr)
library(magrittr)
library(simrabid)
library(sf)
library(raster)
library(lubridate)
library(ggplot2)
library(patchwork)

# study area
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load("data/sd_census_data.rda")
load("data/sd_vacc_data.rda")
load(fp("analysis/out/incursions.csv"))
load("data/sd_case_data.rda")

# source
source("R/sd_data.R")
source("R/utils-data.R")
source("R/figure_funs.R")

# Figure 1 A = time series of cases
sd_case_data %>%
  mutate(date = floor_date(dmy(symptoms_started), unit = "months")) %>%
  filter(year(date) >= 2002 & year(date) <= 2020) %>%
  ggplot() +
  geom_bar(aes(x = date), fill = "#bd2722") +
  scale_x_date(date_breaks = "24 months", date_labels = "%b %Y") + 
  labs(x = "", y = "Number of dog rabies cases", tag = "A") +
  theme_proj -> case_ts_a

# Figure 1B = spatial location of cases
ggplot(sd_shapefile) +
  geom_sf(color = "white", fill = "grey") + 
  geom_point(data = sd_case_data, aes(x = utm_easting, y = utm_northing), 
             color = "#bd2722", alpha = 0.5, size = 1) +
  labs(tag = "B") + 
  theme_map -> case_map_b

# Figure 1D = vaccination timeseries
sd_vacc_data %>% 
  mutate(date = floor_date(dmy(start_date), unit = "year")) %>%
  group_by(date) %>%
  summarize(dogs_vacc = sum(dogs_vaccinated)) %>%
  ggplot() +
  geom_col(aes(x = date, y = dogs_vacc), fill = "#622b67", alpha = 0.75) +
  scale_x_date(date_breaks = "24 months", date_labels = "%Y") +
  labs(x = "", y = "Number dogs vaccinated", tag = "C") +
  theme_proj -> vacc_ts_c

# census data
ggplot(sd_shapefile) +
  geom_sf(color = "white", fill = "grey") + 
  geom_point(data = sd_census_data, 
             aes(x = utm_easting, y = utm_northing, 
                 color = adult_dogs + pups), 
             alpha = 0.25, size = 0.75) +
  scale_color_distiller(palette = "BuPu", trans = "sqrt", direction = 1) +
  labs(tag = "D", color = "Dogs in HH") + 
  theme_map + 
  theme(legend.position="bottom") -> census_data_d

fig1 <- case_ts_a + case_map_b + vacc_ts_c + census_data_d
  
ggsave("analysis/figs/fig_data.jpeg", height = 8, width = 8)
