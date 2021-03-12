
get_observed_data <- function(sd_case_data, mod_specs, mod_start) {
  # unit of simulation
  units <- case_when(mod_specs$days_in_step == 7 ~ "weeks", 
                     mod_specs$days_in_step == 1 ~ "days", 
                     mod_specs$days_in_step == 30.5 ~ "months")
  tmult <- case_when(mod_specs$days_in_step == 7 ~ 52, 
                     mod_specs$days_in_step == 1 ~ 365, 
                     mod_specs$days_in_step == 30.5 ~ 12)
  sd_case_data %<>%
    mutate(t_infectious = get_timestep(symptoms_started, 
                                origin_date = "01-01-2002",
                                date_fun = lubridate::dmy,
                                units = units), 
           month = get_timestep(symptoms_started, 
                                origin_date = "01-01-2002",
                                date_fun = lubridate::dmy,
                                units = "months"),
           cell_id = cellFromXY(mod_start$rast, cbind(utm_easting, utm_northing)), 
           infected = TRUE, detected = TRUE) 
  
  # cases by month & cell
  cases_by_month <- tabulate(sd_case_data$month, nbins = mod_specs$nyears * 12)
  cases_by_cell <- tabulate(sd_case_data$cell_id, nbins = ncell(mod_start$rast))
  
  # set-up for input into summ_stats
  ncells <- ncell(mod_start$rast)
  tmax <- mod_specs$nyears * tmult
  extra_pars <- list(obs_data =  list(cases_by_month = cases_by_month,
                                      cases_by_cell = cases_by_cell))
  days_in_step <- mod_specs$days_in_step
  
  sd_case_data %>%
    dplyr::select(t_infectious, x_coord = utm_easting, y_coord = utm_northing, 
           infected, detected, cell_id) %>%
    as.data.table() -> I_dt
  
  # apply summary functions to the observed data
  obs_sstats <- inc_stats()
  
  return(list(obs_data = extra_pars$obs_data, obs_sstats = obs_sstats))
}
