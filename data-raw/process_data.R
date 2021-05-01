# ------------------------------------------------------------------------------
#' Proccess data
# ------------------------------------------------------------------------------

library(readr)
library(dplyr)
library(sf)
library(raster)
library(fasterize)
library(data.table)
library(janitor)
library(magrittr)
library(lubridate)
select <- dplyr::select
source("R/utils-data.R")

# Data -------------------------------------

# load in shapefile & other data
sd_shapefile <- st_read(system.file("extdata/sd_shapefile.shp", 
                                    package = "simrabid"))
# district wide census (@ vill level)
sd_census <- read_csv(get_latest("data-raw/wisemonkey", "Census"))
case_dt <- read_csv(get_latest("data-raw/wisemonkey", "Animal_Contact_Tracing"))
vacc_dt <- read_csv(get_latest("data-raw/wisemonkey", "Vaccination"))
inc_dt <- read_csv("data-raw/incursions.csv")
sd_pops <- read_csv("data-raw/SerengetiPop.csv")

# Clean census data and out neccessary bits -----
sd_census %>%
  clean_names() %>%
  select(village_2002, utm_easting, utm_northing,
         adults, children, adult_dogs = dogs_3_months,
         adult_dogs_vacc = vaccinated_dogs_3_months,
         pups = pups_3_months, pups_vacc = vaccinated_pups_3_months,
         date = actvitity_date) %>%
  filter(!is.na(utm_easting), !is.na(utm_northing)) -> sd_census_data

usethis::use_data(sd_census_data, overwrite = TRUE)

# Case data ----
bounds <- st_bbox(sd_shapefile)

# First clean data for tree building (i.e. to identify incursions)
case_dt %>%
  clean_names() %>%
  select(id, biter_id, 
         village_2002, utm_easting, utm_northing, biter_id, species, suspect,
         symptoms_started_known, symptoms_started,
         symptoms_started_accuracy)  %>%
  filter(suspect %in% "Yes", symptoms_started_known,
         species %in% "Domestic dog", 
         year(dmy(symptoms_started)) >= 2002, 
         year(dmy(symptoms_started)) <= 2020, 
         utm_easting >= bounds$xmin & utm_easting <= bounds$xmax, 
         utm_northing >= bounds$ymin & utm_northing <= bounds$ymax, 
         !is.na(id)) -> sd_case_data_trees # na id's are ones that do not fall within the district

usethis::use_data(sd_case_data_trees, overwrite = TRUE)

# Exclude duplicate ids for case data (for spatial stats & time series)
sd_case_data_trees %>%
  group_by(id) %>%
  slice(1) %>% # filter out duplicate ids
  ungroup() -> sd_case_data

usethis::use_data(sd_case_data, overwrite = TRUE)

# Vacc data ----
vacc_dt %>%
  clean_names() %>%
  select(village_2002, start_date, end_date,
         doses, dogs_vaccinated) %>%
  # when dogs_vaccinated not specified assume 95% of doses go to dogs
  mutate(dogs_vaccinated = coalesce(dogs_vaccinated,
                                         round(doses * 0.95))) -> sd_vacc_data

correct_names <- tribble(~village_2002, ~corrected,
                          "Mbirikili", "Bonchugu",
                          "Kebanchabancha", "Kebanchebanche",
                          "Natta", "Mbisso",
                          "Natta Mbiso", "Mbisso",
                          "Nyamisingisi", "Nyamasingisi")

sd_vacc_data %<>%
  left_join(correct_names) %>%
  mutate(village_2002 = coalesce(corrected, village_2002),
        match_2002 = sd_pops$VILLCODES[match(village_2002,
                                              sd_pops$Village_2002)],
        match_2012 = sd_pops$VILLCODES[match(village_2002,
                                               sd_pops$Village_2012)],
        villcode = coalesce(match_2002, match_2012)) %>%
  select(-contains(c("match", "corrected")))

usethis::use_data(sd_vacc_data, overwrite = TRUE)

# incursions -----
inc_dt %>%
  clean_names() %>%
  filter(incursions == TRUE) %>%
  select(id, date = d_symptoms) %>%
  left_join(select(clean_names(case_dt), id, x_coord = utm_easting, 
                   y_coord = utm_northing)) -> incursions
usethis::use_data(incursions, overwrite = TRUE)
