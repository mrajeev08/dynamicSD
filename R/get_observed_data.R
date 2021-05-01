
get_observed_data <- function(sd_case_data, cand, out) {
  
  days_in_step <- cand$days_in_step
  start_date <- cand$start_date
  
  # get timesteps & tmax using approximate dates
  step <- 365.25 / days_in_step
  tmax <- ceiling(as.numeric(lubridate::ymd(cand$apprx_end_date) - lubridate::ymd(start_date)) / days_in_step)
  
  # actual end date (to ensure same number of days in each timestep)
  end_date <- lubridate::ymd(start_date) + tmax * days_in_step
  end_month <- get_timestep(end_date, origin_date = start_date, 
                            date_fun = lubridate::ymd, 
                            days_in_step = 30.5)
  sd_case_data %<>%
    mutate(t_infectious = get_timestep(symptoms_started, 
                                       origin_date = start_date, 
                                       date_fun = lubridate::dmy, 
                                       days_in_step = days_in_step), 
           month = get_timestep(symptoms_started, 
                                origin_date = start_date,
                                date_fun = lubridate::dmy,
                                days_in_step = 30.5),
           cell_id = cellFromXY(out$rast, cbind(utm_easting, utm_northing)), 
           infected = TRUE, detected = TRUE) 
  
  # cases by month & cell
  cases_by_month <- tabulate(sd_case_data$month, nbins = end_month)
  cases_by_cell <- tabulate(sd_case_data$cell_id, nbins = ncell(out$rast))

  # set-up for input into summ_stats
  ncells <- ncell(out$rast)
  extra_pars <- list(obs_data =  list(cases_by_month = cases_by_month,
                                      cases_by_cell = cases_by_cell))
  t <- tmax
  prop_start_pop <- 7e4/sum(out$start_pop, na.rm = TRUE)
  break_threshold <- 0
  
  sd_case_data %>%
    dplyr::select(t_infectious, x_coord = utm_easting, y_coord = utm_northing, 
                  infected, detected, cell_id) %>%
    as.data.table() -> I_dt
  
  # apply summary functions to the observed data (will also generate 5 noise stats)
  obs_sstats <- inc_stats()
  
  return(list(obs_data = extra_pars$obs_data, obs_sstats = obs_sstats))
}
