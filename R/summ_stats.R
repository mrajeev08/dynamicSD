# Summary stat functions
# To try
# max # of cases
# mean # of cases
# median # of cases
# pop growth overall (at district scale) (so that its okay if timed out)
# acf (1 - 10)
# mean distance between cases at 1 month (Rcpp)
# this mean normalized by mean overall (Rcpp)
# spatial loss function (mismatch between where cases are?)
# ks test statistic for frequency of case counts per month (temporal bit) (so that its okay if timed out)


# I_dt <- data.table(tstep = runif(1e4, 0, 4*24),
#                    x_coord = runif(1e4, 0, 1),
#                    y_coord = runif(1e4, 0, 1))
#

# Incidence stats
# data = inc_hist + breaks (which should be 1 longer than length of inc_hist)
# data = # of cases per cell (order by the cell_id)
# do this with example of different sizes & benchmark
inc_stats <- function(names = c("I_dt", "ncells", "tmax", "extra_pars", 
                                "days_in_step")) {

  # Get the objects you need from the environment above this one
  list2env(use_mget(names, envir_num = 2), envir = environment())

  # filter I_dt to successful transmission events & detected cases
  I_dt <- I_dt[infected & detected]
  
    
  # aggregate cols by timestep
  ncols_sum <- floor(30.5 / days_in_step)
  
  # Summarize monthly cases
  I_ts <- tabulate(floor(I_dt$t_infectious), tmax)
  I_ts <- sum_to_month(I_total, nc = ncols_sum)
  
  # return stats on the time series (should not have NAs)
  max_I <- max(I_ts)
  median_I <- median(I_ts)
  mean_I <- mean(I_ts)
  
  # temporal corr
  acfs <- as.vector(acf(I_ts, lag.max = 10, plot = FALSE)$acf)[-1]
  names(acfs) <- paste0("acf_lag", 1:10)
  
  # spatial corr
  normalized <- mean(dist(cbind(I_dt$x_coord, I_dt$y_coord)))
  mean_dist_4wks <- get_mean_dist(t_dt = I_dt[, .(tstep, x_coord, y_coord)],
                                  t_window = 8, samp_max = 1e4)
  mean_dist_4wks_norm <- mean_dist/normalized
  mean_dist_8wks <- get_mean_dist(t_dt = I_dt[, .(tstep, x_coord, y_coord)],
                                  t_window = 8, samp_max = 1e4)
  mean_dist_8wks_norm <- mean_dist/normalized
  
  # temporal loss
  temp_rmse <- sqrt(mean((I_ts - data$cases_by_month)^2))
  temp_ss <- sum((I_ts - data$cases_by_month)^2)
  
  # spatial loss
  I_cell <- tabulate(I_dt$cell_id, ncells)
  spat_rmse <- sqrt(mean((I_cell - data$cases_by_cell)^2))
  spat_ss <- sum((I_cell - data$cases_by_cell)^2)
  spat_loss <- mean(abs(I_cell - data$cases_by_cell))
  
  # ks discrete statistic
  inc_hist <- hist(I_ts, breaks = data$breaks)$count
  ks_stat <- max(abs(inc_hist - data$inc_hist))
  hist_ss <- sum((I_cell/sum(I_cell) - data$cases_by_cell/sum(data$cases_by_cell))^2)
  hist_rmse <- sqrt(mean((I_cell/sum(I_cell) - data$cases_by_cell/sum(data$cases_by_cell))^2))
  
  
  return(c(list(max_I = max_I, median_I = median_I, mean_I = mean_I,
                ks_stat = ks_stat,
                hist_ss = hist_ss, hist_rmse = hist_rmse,
                mean_dist_4wks = mean_dist_4wks,
                mean_dist_8wks = mean_dist_8wks,
                mean_dist_4wks_norm = mean_dist_4wks_norm,
                mean_dist_8wks_norm = mean_dist_8wks_norm, 
                spat_rmse = spat_rmse,
                spat_ss = spat_ss, 
                temp_rmse = temp_rmse,
                temp_ss = temp_ss), as.list(acfs)))

}


# Fast acf

# Mean distance between cases within one month of each other
# (lag in 1 direction only to avoid mutliple reps?)

# Normalized by mean distance between all cases (dist_mat mean)

# Naive version with coords (beyond 10,000 coords this becomes very slow)
# tstep <- runif(1e4, 0, 4*24)
# x_coord <- runif(1e4, 0, 1)
# y_coord <- runif(1e4, 0, 1)
# mean_dist_window(x_coord, y_coord, tstep, 4)

mean_dist_window <- function(I_dt, t_window = 4,
                             normalize = TRUE) {

  mu_dist <- rep(NA, length(x_coord))

  for(i in seq_len(length(x_coord))) {

    diff_t <- tstep[i] - tstep

    # filter to ones with cases in preceding twindow
    within <- diff_t > 0 & diff_t < t_window

    if(sum(within) > 0 ) {
      mu_dist[i] <- mean((x_coord[i] - x_coord[within])^2 + (y_coord[i] - y_coord[within])^2)
    } else {
      next
    }

  }

  mean_dist <- mean(mu_dist, na.rm = TRUE)

  if(normalize) {
    mean_dist <- mean_dist/mean(dist(cbind(x_coord, y_coord)))
  }

  return(mean_dist)

}


mean_dist_dt<- function(I_dt, t_window = 4,
                             normalize = TRUE) {

  mu_dist <- rep(NA, length(x_coord))

  for(i in seq_len(length(x_coord))) {

    diff_t <- tstep[i] - tstep

    # filter to ones with cases in preceding twindow
    within <- diff_t > 0 & diff_t < t_window

    if(sum(within) > 0 ) {
      mu_dist[i] <- mean((x_coord[i] - x_coord[within])^2 + (y_coord[i] - y_coord[within])^2)
    } else {
      next
    }

  }

  mean_dist <- mean(mu_dist, na.rm = TRUE)

  if(normalize) {
    mean_dist <- mean_dist/mean(dist(cbind(x_coord, y_coord)))
  }

  return(mean_dist)

}

# Get the x & ycoords for up to N lags & then filter to the diff in times and take the means
# create a data.table with the tsteps
# this is the fastest!
# set it to do max 10k samples
dist_fun <- function(x_coord, y_coord, x_coord_to, y_coord_to) {

 mean(sqrt((x_coord - x_coord_to)^2 + (y_coord - y_coord_to)^2))

}

get_mean_dist <- function(t_dt = I_dt[, .(tstep, x_coord, y_coord)],
                          t_window = 4, samp_max = 1e4) {

  if(nrow(t_dt) > samp_max) {
    t_dt <- t_dt[sample(.N, samp_max)]
  }

  c_dt <- t_dt[ , .(max = tstep,
                      min = tstep - 4,
                      x_coord_to = x_coord,
                      y_coord_to = y_coord)]

  new <- t_dt[c_dt,
                    on = .(tstep < max, tstep >= min),
                    allow.cartesian = TRUE, nomatch = NULL]

  return(dist_fun(new$x_coord, new$y_coord, new$x_coord_to, new$y_coord_to))

}

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

