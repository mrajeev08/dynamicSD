# Time series from simulations
ts_stats <- function(names = c("I_dt", "extra_pars", 
                               "start_date", "days_in_step", 
                               "S_mat", "N_mat",
                               "prop_start_pop", 
                               "break_threshold")) {
  
  # Get the objects you need from the environment above this one
  list2env(use_mget(names, envir_num = 2), envir = environment())
  
  obs_data <- extra_pars$obs_data
  
  # filter I_dt to successful transmission events & detected cases
  I_dt <- I_dt[infected & detected]
  
  # aggregate cols by timestep (dates!)
  I_dt[, date := get_date(origin_date = start_date, 
                          tstep = t_infectious, 
                          days_in_step = days_in_step)]
  I_dt[, cal_month := get_cal_month(date,
                                    origin_date = start_date)]
  
  # Summarize monthly cases
  cal_month_max <- length(obs_data$cases_by_month)
  I_ts <- tabulate(I_dt$cal_month, cal_month_max)
  
  # Summarize S/N monthly (@ district level)
  ts_dates <- get_date(origin_date = start_date, 
                       tstep = 1:ncol(S_mat), 
                       days_in_step = days_in_step)
  ts_months <- get_cal_month(ts_dates, origin_date = start_date)
  ts_select <- match(seq_len(cal_month_max), ts_months)
  S <- colSums(S_mat)[ts_select]
  N <- colSums(N_mat)[ts_select]
  cov <- S/N # rough estimate
  
  # Summarize monthly introduced cases
  incs_ts <- tabulate(I_dt$cal_month[I_dt$progen_id == -1], cal_month_max)
  incs_success <- tabulate(
    I_dt$cal_month[I_dt$progen_id == -1 & I_dt$id %in% I_dt$progen_id], 
    cal_month_max)
  
  stopped <- prop_start_pop < break_threshold
  
  # out data.table
  return(data.table(cal_month = seq_len(cal_month_max),
                    I_ts = I_ts, incs_success = incs_success, 
                    incs_ts = incs_ts, cov = cov, S = S, N = N,
                    stopped = stopped))
  
}

# Summary stats for abc ----
inc_stats <- function(names = c("I_dt", "ncells", "tmax", "extra_pars", 
                                "days_in_step", "prop_start_pop", 
                                "break_threshold", "t", "start_date")) {

  # Get the objects you need from the environment above this one
  list2env(use_mget(names, envir_num = 2), envir = environment())
  obs_data <- extra_pars$obs_data
  
  # filter I_dt to successful transmission events & detected cases
  I_dt <- I_dt[infected & detected]
  
  # aggregate cols by timestep (dates!)
  I_dt[, date := get_date(origin_date = start_date, 
                          tstep = t_infectious, 
                          days_in_step = days_in_step)]
  
  I_dt[, ts_agg := get_timestep(date, 
                                origin_date = start_date,
                                date_fun = lubridate::ymd,
                                days_in_step = 30.5)]
  I_dt[, ts_agg_week := get_timestep(date, 
                                     origin_date = start_date,
                                     date_fun = lubridate::ymd,
                                     days_in_step = 7)]
  
  # Weekly acfs
  I_ts_week <- tabulate(I_dt$ts_agg_week, max(I_dt$ts_agg_week))
  # temporal corr
  acfs_week <- as.vector(acf(I_ts_week, lag.max = 10, plot = FALSE,
                             na.action = na.pass)$acf)[-1]
  
  if(length(acfs_week) < 10) {
    acfs_week <- c(acfs_week, rep(0, 10 - length(acfs_week)))
  }
  names(acfs_week) <- paste0("acf_week", 1:10)
  
  # Summarize monthly cases
  I_ts <- tabulate(I_dt$ts_agg, length(obs_data$cases_by_month))
  
  # account for simulations stopped early
  stopped <- prop_start_pop < break_threshold
  
  if(stopped) {
    brk_date <- get_date(origin_date = start_date, 
                         tstep = t, 
                         days_in_step = days_in_step)
    brk <- get_timestep(date = brk_date, 
                        origin_date = start_date,
                        date_fun = lubridate::ymd,
                        days_in_step = 30.5)
    I_ts <- I_ts[1:brk]
  }
  
  # return stats on the time series (should not have NAs)
  max_I <- max(I_ts)
  median_I <- median(I_ts)
  mean_I <- mean(I_ts)
  
  # temporal corr
  acfs <- as.vector(acf(I_ts, lag.max = 6, plot = FALSE,
                        na.action = na.pass)$acf)[-1]
  if(length(acfs) < 6) {
    acfs <- c(acfs, rep(0, 6 - length(acfs)))
  }
  names(acfs) <- paste0("acf_month", 1:6)
  
  # spatial corr
  I_dt <- I_dt[!is.na(x_coord) | !is.na(y_coord)]
  mean_dist_all <- mean(dist(cbind(I_dt$x_coord, I_dt$y_coord)))
  mean_dist_4wks <- get_mean_dist(t_dt = I_dt[, .(t_infectious, x_coord, y_coord)],
                                  t_window = 4, samp_max = 1e4)
  mean_dist_4wks_norm <- mean_dist_4wks/mean_dist_all
  mean_dist_8wks <- get_mean_dist(t_dt = I_dt[, .(t_infectious, x_coord, y_coord)],
                                  t_window = 8, samp_max = 1e4)
  mean_dist_8wks_norm <- mean_dist_8wks/mean_dist_all
  
  # temporal loss
  if(stopped) {
    I_ts <- sample(I_ts, length(obs_data$cases_by_month), replace = TRUE)
  }
  
  temp_rmse <- sqrt(mean((I_ts - obs_data$cases_by_month)^2))
  temp_loss <- mean(abs(I_ts - obs_data$cases_by_month))
  
  # spatial loss
  I_cell <- tabulate(I_dt$cell_id, ncells)
  obs_data_prop <- obs_data$cases_by_cell/sum(obs_data$cases_by_cell)
  I_cell_prop <- I_cell/sum(I_cell)
  spat_rmse <- sqrt(mean((I_cell_prop -  obs_data_prop)^2))
  spat_loss <- sum(abs(I_cell_prop - obs_data_prop))

  # ks discrete statistic
  kb_temp <- kb_stat(obs_data$cases_by_month, I_ts)
  kb_spat <- kb_stat(obs_data$cases_by_cell, I_cell)
  
  # out data.table
  return(data.table(max_I, median_I, mean_I, kb_temp, kb_spat,
                    mean_dist_4wks, mean_dist_8wks, 
                    mean_dist_4wks_norm, mean_dist_8wks_norm, 
                    spat_rmse, spat_loss,
                    temp_rmse, temp_loss, 
                    mean_dist_all,
                    t(acfs), t(acfs_week), 
                    prop_start_pop, break_threshold, 
                    stopped)) 

}


# Get mean distance between cases N week lag before index case
get_mean_dist <- function(t_dt = I_dt[, .(t_infectious, x_coord, y_coord)],
                          t_window = 4, samp_max = 1e4) {

  if(nrow(t_dt) > samp_max) {
    t_dt <- t_dt[sample(.N, samp_max)]
  }

  c_dt <- t_dt[ , .(max = t_infectious,
                      min = t_infectious - t_window,
                      x_coord_to = x_coord,
                      y_coord_to = y_coord)]

  new <- t_dt[c_dt,
                    on = .(t_infectious < max, t_infectious >= min),
                    allow.cartesian = TRUE, nomatch = NULL]

  return(dist_fun(new$x_coord, new$y_coord, new$x_coord_to, new$y_coord_to))

}

# helper for distance fun
dist_fun <- function(x_coord, y_coord, x_coord_to, y_coord_to) {
  
  mean(sqrt((x_coord - x_coord_to)^2 + (y_coord - y_coord_to)^2))
  
}

# Function for benchmarking + testing simulation output
# Output = 1 row per sim
test_sim <- function(names = c("I_dt", "S_mat", "N_mat", "tmax", "days_in_step")) {
  
  start <- Sys.time()
  
  # Get the objects you need from the environment above this one
  list2env(use_mget(names, envir_num = 2), envir = environment())
  
  # Filter to infected (this also means it will no longer point to I_dt)
  I_dt <- I_dt[infected == TRUE]
  
  # aggregate cols by timestep
  ncols_sum <- floor(30.5 / days_in_step)
  if(days_in_step == 1) denom <- 7 else denom <- 1
  
  I_total <- tabulate(floor(I_dt$t_infectious), tmax)
  I_detected <- tabulate(floor(I_dt$t_infectious[I_dt$detected == TRUE]), tmax)
  I_local <- tabulate(floor(I_dt$t_infectious[I_dt$progen_id > 0]), tmax)
  
  # Summarize monthly cases
  I_total <- sum_to_month(I_total, nc = ncols_sum)
  I_local <- sum_to_month(I_local, nc = ncols_sum)
  I_detected <- sum_to_month(I_detected, nc = ncols_sum)
  
  # Summarize S/N monthly
  S <- inds_to_month(colSums(S_mat), nc = ncols_sum)
  N <- inds_to_month(colSums(N_mat), nc = ncols_sum)
  cov <- S/N
  
  # distances
  if(nrow(I_dt) > 0) {
    
    mean_dist_m <- mean(dist_linked(I_dt))
    max_dist_m <- max(dist_linked(I_dt))
    
    # times
    mean_times_wks <- mean(times_linked(I_dt))/denom
    max_times_wks <- max(times_linked(I_dt))/denom
    
  } else {
    
    mean_dist_m <- max_dist_m <- mean_times_wks <- max_times_wks <- NA
    
  }
  
  
  # data.table with tstep as additional covariate
  data.table(max_monthinc_total = max(I_total),
             max_monthinc_local = max(I_local), 
             max_monthinc_detected = max(I_detected),
             max_sus = max(S), 
             max_N = max(N), 
             max_cov = max(cov), 
             mean_dist_m, 
             max_dist_m, 
             mean_times_wks, 
             max_times_wks, 
             summ_time = as.numeric(Sys.time() - start))
}

# summarize every 4 columns
sum_to_month <- function(ts, nc = 4) {
  
  if(!is.null(dim(ts))) {
    ts <- colSums(ts)
  }
  nv <- length(ts)
  if (nv %% nc)
    ts[ceiling(nv / nc) * nc] <- NA
  colSums(matrix(ts, nc), na.rm = TRUE)
  
}

# summarize every 4 columns
inds_to_month <- function(ts, nc = 4) {
  
  maxl <- floor(length(ts) / nc)
  ts <- ts[seq(1, length(ts), by = nc)]
  ts[1:maxl]
  
}

# average euclidean distance between linked cases
dist_linked <- function(I_dt) {
  
  coords <- I_dt[, c("x_coord", "y_coord", "id", "progen_id")]
  coords <- coords[coords, on = c("id" = "progen_id")][!is.na(progen_id)]
  coords[, dist_m := sqrt((x_coord - i.x_coord)^2 + (y_coord - i.y_coord)^2)]
  return(coords$dist_m)
}

# average times between linked cases
times_linked <- function(I_dt) {
  
  times <- I_dt[t_infected > 0][, ts := t_infectious - t_infected]
  
  return(times$ts)
}

kb_stat <- function(x, y) {
  
  pmax <- max(c(x, y))
  px <- tabulate(x, pmax)/sum(x)
  py <- tabulate(y, pmax)/sum(y)
  
  sum(exp(px) * (log(exp(px) / exp(py))))

}
