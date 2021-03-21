run_simrabid <- function(cand, 
                         mod_specs,
                         param_ests,
                         param_defaults,
                         nsims, 
                         extra_pars,
                         vacc_dt,
                         combine_fun, 
                         summary_fun, 
                         merge_fun = list,
                         secondary_fun,
                         weight_covars, 
                         weight_params, 
                         multi = FALSE, 
                         convolve_steps = TRUE, 
                         sim_vacc = "none", 
                         ...) {
  
  
  start_up <- setup_sim(start_date = cand$start_date, 
                        apprx_end_date =  cand$apprx_end_date, 
                        days_in_step = cand$days_in_step, 
                        rast = mod_specs$rast,
                        death_rate_annual = mod_specs$death_rate_annual,
                        birth_rate_annual = mod_specs$birth_rate_annual,
                        waning_rate_annual = cand$waning,
                        block_fun = block_cells,
                        params = list(start_pop = mod_specs$start_pop),
                        by_admin = cand$by_admin)
  
  if(!cand$estincs) {
    
    # add cell ids & timesteps of incursion to parameter defaults
    incursions %<>%
      mutate(cell_id = raster::cellFromXY(mod_specs$rast,
                                          cbind(x_coord, y_coord)), 
             tstep = get_timestep(date, 
                                  origin_date = cand$start_date,
                                  date_fun = lubridate::ymd, 
                                  days_in_step = cand$days_in_step))
    
    param_defaults <- c(param_defaults, 
                        list(cell_ids = incursions$cell_id, 
                             tstep = incursions$tstep))
  }
  
  if(cand$weights) {
    weights <- cell_weights(covars = weight_covars,
                            params = weight_params,
                            start_up$ncell,
                            leave_bounds = cand$leave_bounds,
                            allow_invalid = cand$allow_invalid,
                            cells_block = start_up$cells_block,
                            cells_out_bounds = start_up$cells_out_bounds)
  } else {
    weights <- NULL
  }
  
  # decide parameters based on args from cand
  disp_fn <- ifelse(cand$sequential, 
                    ifelse(convolve_steps, 
                           function(n, params) { # actually steps are twice as long on avg
                             rweibull(n, shape = params$steps_shape, 
                                      scale = params$steps_scale) + 
                              rweibull(n, shape = params$steps_shape, 
                                       scale = params$steps_scale)
                           }, simrabid::steps_weibull), 
                    simrabid::dispersal_lognorm)
  
  inc_fn <- ifelse(cand$estincs, sim_incursions_pois, sim_incursions_hardwired)
  move_fn <- ifelse(cand$weights, simrabid::sim_movement_prob, simrabid::sim_movement_continuous)
  
  # removing args if they're included in the priors
  param_defaults <- param_defaults[!(names(param_defaults) %in% names(param_ests))]
  
  # check priors & make sure they're reproducible
  set.seed(cand$seed)
  
  param_ests <- lapply(param_ests, function(x) {
    if(is.function(x)) {
      x(nsims)
    } else {
      if (length(x) == nsims) {
        x
      } else {
        if(length(x) == 1) {
          rep(x, nsims)
        } else {
          stop("Parameter passed, but it is not a function, of length 1, or of length nsims!")
        }
      }
    }
  })
  
  # functions to export to foreach loop
  export_funs <- c(list_funs("R/utils-data.R"), 
                   list_funs("R/summ_stats.R"), 
                   list_funs("R/conn_metrics.R")) # this should be an argument
  
  out_sims <-
    foreach(j = seq_len(nsims), 
            .combine = combine_fun, 
            .multicombine = multi, 
            .packages = c("simrabid", # this should be an argument
                          "data.table", 
                          "raster", 
                          "magrittr", 
                          "igraph"),
            .export = export_funs,
            .options.RNG = cand$seed) %dorng% {
              
              # draw a val 
              pars <- lapply(param_ests, function(x) x[j])
              
              out <- 
                tryCatch(
                  expr = {
                    
                    if(sim_vacc == "fixed") {
                      
                      # same vills for each sim
                      max_id <- max(start_up$loc_ids)
                      set.seed(cand$seed * j)
                      locs <- seq_len(max_id)[rbinom(max_id, size = 1, prob = vacc_dt$vacc_prop)]
                      vacc_dt_i <- sim_campaigns(locs = locs,
                                               campaign_prob = 1, 
                                               coverage = vacc_dt$vacc_cov, 
                                               sim_years = vacc_dt$years, 
                                               burn_in_years = vacc_dt$burn_in) 
                      cover <- TRUE
                    }
                    
                    if(sim_vacc == "random") {
                      
                      vacc_dt_i <- sim_campaigns(locs = seq_len(max(start_up$loc_ids)),
                                               campaign_prob = vacc_dt$vacc_prop, 
                                               coverage = vacc_dt$vacc_cov, 
                                               sim_years = vacc_dt$years, 
                                               burn_in_years = vacc_dt$burn_in)
                      cover <- TRUE
                      
                    }
                    
                    if(sim_vacc == "none") {
                      cover <- FALSE
                      vacc_dt_i <- vacc_dt
                    }
                    
                    simstats <- simrabid(start_up, 
                                         start_vacc = cand$start_vacc, 
                                         I_seeds = cand$I_seeds, 
                                         vacc_dt = vacc_dt_i,
                                         params = c(pars, 
                                                    param_defaults),
                                         days_in_step = cand$days_in_step,
                                         observe_fun = beta_detect_monthly,
                                         serial_fun = serial_lognorm,
                                         dispersal_fun = disp_fn, 
                                         secondary_fun = secondary_fun, # function argument 
                                         incursion_fun = inc_fn, 
                                         movement_fun = move_fn,
                                         sequential = cand$sequential, 
                                         allow_invalid = cand$allow_invalid,
                                         leave_bounds = cand$leave_bounds, 
                                         max_tries = 100,
                                         summary_fun = summary_fun, # function argument
                                         track = cand$track,
                                         weights = weights, 
                                         row_probs = NULL,
                                         coverage = cover, # this should be an argument!
                                         break_threshold = cand$break_threshold,
                                         by_admin = cand$by_admin,
                                         extra_pars = extra_pars, 
                                         ...)
                    
                  
                    out <- merge_fun(simstats, t(pars), sim = j)
                  },
                  error = function(e){
                    # write out the error & the parameter log
                    err <- data.table(t(pars), cand,  
                                      message = e$message)
                    fwrite(err, 
                              paste0("logs/run_simrabid_err", 
                                     format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                     ".csv"),
                              row.names = FALSE)
                    NULL
                  }
                )
                
            }
  
  if(!all(out_sims$sim %in% 1:nsims)) {
    message("Warning: some simulations failed, check the logs directory for
            more info!")
  }
  
  if(is.null(out_sims)) {
    message("No simulations were successful! Check the logs for more info, 
            returning NA for now.")
    out_sims <- list(NA)
  }
  
  return(out_sims)
  
}
