# stats for paper ----

library(here)
library(data.table)
library(dplyr)
library(magrittr)
library(simrabid)
library(sf)
library(raster)
library(lubridate)

# study area
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
load(here("data/sd_census_data.rda"))
load(here("data/sd_vacc_data.rda"))
load(here("data/incursions.rda"))
load(here("data/sd_case_data.rda"))

# source
source(here("R/sd_data.R"))
source(here("R/utils-data.R"))

# stats on human pop and the district
human_pop_2012 <- round(sum(sd_shapefile$pop_2012), -3)
sd_size_km2 <- round(sqrt(as.numeric(st_area(st_combine(sd_shapefile)))/1000), -2)

# get the startup space
out <- get_sd_pops(sd_shapefile, res_m = 1000,
                   sd_census_data, death_rate_annual = 0.465)

# Group by villcode & get human to dog ratios by village
sd_census_data %<>%
  mutate(villcode = get_ids(x_coord = utm_easting,
                            y_coord = utm_northing,
                            rast = out$rast,
                            id_col = sd_shapefile$villcode))

sd_census_data %>%
  group_by(villcode) %>%
  tidyr::replace_na(list(adults = 0, children = 0,
                         adult_dogs = 0, pups = 0)) %>%
  summarize(human_pop = sum(adults + children, na.rm = TRUE),
            dog_pop = sum(adult_dogs + pups, na.rm = TRUE)) %>%
  mutate(hdr = human_pop/dog_pop) -> sd_hdr
sd_hdr_range <- paste(ceiling(range(sd_hdr$hdr)), collapse = " to ")

# out stats -----
stats <- list(human_pop_2012 = human_pop_2012, 
              sd_hdr_range = sd_hdr_range, 
              sd_size_km2 = sd_size_km2, 
              dates_census = paste0(year(dmy(range(sd_census_data$date))), collapse = " - "))

              