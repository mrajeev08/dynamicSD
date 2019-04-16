## You need the 
sim.IBM <- function(grid = as.matrix(SD_raster), data = grid_data, vacc = vacc_mat,
                    I_seed = 2, start_vacc = 0.2, R_0 = 1.1, k = 0.2, 
                    inc_rate = incs, # weekly number of incursions and scaling factor 
                    sigma = get.prob(rate = 7/22.3, step = 1), # weekly rate to prob exp -> inf
                    births = get.prob(rate = 0.44 + (grid_data$growth - 1), step = 52), # annual birth rate to prob
                    cell_id = grid_data$cell_id, 
                    mu = get.prob(rate = 0.44, step = 52), # annual death rate to prob 
                    ntimes = tmax, # maximum time step,
                    nu = get.prob(rate = 0.33, step = 52), # annual waning rate to prob,
                    p_revacc = 0.5, # probability of revaccination
                    p_vacc = 1, # probability of vacc by grid cell
                    dispersalShape = 0.3484, dispersalScale = 41.28/100, # from Rebecca
                    x_topl = bbox(SD_raster)[1, "min"], 
                    y_topl = bbox(SD_raster)[2, "max"], res_m = res(SD_raster)[1], 
                    return_coords = FALSE, vacc_sim = "empirical", cov_val = 1, 
                    perc_val = 1,...) {
  
  # 1. Set-up and tstep 1 -----------------------------------------------------------------------
  nlocs = nrow(data)
  row_id <- 1:nlocs
  
  # state matrices
  S <- V <- N <- matrix(0, nrow = nlocs, ncol = ntimes)
  I_dist <- I_all <- E <- matrix(0, nrow = nlocs, ncol = ntimes)
  
  ## Starting pop + sus 
  N[, 1] <- data$start_pop
  V[, 1] <- rbinom(n = nlocs, size = data$start_pop, prob = start_vacc)
  S[, 1] <- N[, 1] - V[, 1]
  
  # 2. Simulating infection in first timestep ------------------------------------------------------
  ##' Populated cells at beginning of simulation 
  ##' (if you want to limit where doggies go in the biting steps, then need to put this within the 
  ##' time loop and have an escape clause for the while loop when landscape is sparse, otherwise
  ##' will spend loads of time moving around trying to find a spot with a doggo!)
  cells_pop <- cell_id[which(N[, 1] > 0)] ## populated cells at beginning
  rows_pop <- row_id[which(N[, 1] > 0)]
  
  # Seeding cases (as if from outside!)
  I_locs <- sample((row_id)[which(!is.na(data$start_pop))], I_seed)   # pick random location to start I seeds
  I_dist[I_locs, 1] <- 1
  
  ## Create infected coordinate data frame
  I_coords <- data.table(ID = 1:I_seed, tstep = rep(1, I_seed), x_coord = data$x_coord[I_locs], 
                         y_coord = data$y_coord[I_locs], progen_ID = rep(0, I_seed),
                         path_ID = rep(NA, I_seed), cell_id = cell_id[I_locs], 
                         trans = NA_integer_, infectious = rep(1, I_seed), 
                         secondaries = 1) ## seeding as single case in first timestep
  
  ## Exposed in first timestep
  counter <- max(I_coords$ID) ## for keeping track of ids
  E_coords_first <- sim.bites(secondaries = I_coords$secondaries, x_coord = I_coords$x_coord, 
                        y_coord = I_coords$y_coord, ids = I_coords$ID,
                        grid = grid, cells_pop = cells_pop, counter = counter, 
                        dispersalShape = dispersalShape, dispersalScale = dispersalScale, 
                        x_topl = x_topl, y_topl = y_topl, res_m = res_m, tstep = 1)
  counter <- E_coords_first[["counter"]]
  E_coords <- E_coords_first[["E_coords_now"]]
  
  if(nrow(E_coords > 0)) {
    E_coords <- sim.trans(E_coords, cell_id = cell_id, row_id = row_id, S = S[, 1], N = N[, 1]) 
    if (nrow(E_coords) > 0) {
      E_count <- E_coords[, list(count = length(ID)), by = cell_id]
      E[row_id[match(E_count$cell_id, cell_id)], 1] <- E_count$count ## add in newly exposed
    }
  }
  
  ## Finish subtracting everyone out
  S[, 1] <- S[, 1] - E[, 1] ## subtract out newly exposed guys from S
  
  ## Timesteps to simulate vaccination for
  vstep <- which(colSums(vacc) > 0, arr.ind = TRUE)

  # 3. Simulate for rest of time steps -------------------------------------------------------------
  for (t in 2:ntimes) {
    # t = 2
    # print(paste(t, "/", tmax, "weeks"))
    
    # 3a. Simulate vaccination ---------------------------------------------------------------------
    ## V transitions, sequetially for competing probs (-waning, - deaths)
    V[, t] <- V[, t-1] - rbinom(nlocs, size = V[, t - 1], prob = mu) ## die first
    waning <- rbinom(nlocs, size = V[, t], prob = nu) ## those that don't die then can be lost to waning
    V[, t] <- V[, t] - waning
    
    ## Susceptible transitions die and also adding in waning vaccinated (so can be revacced)
    S[, t] <- S[, t-1] - rbinom(nlocs, size = S[, t-1], prob = mu)  ## die first
    ## then add in waning + new borns
    S[, t] <- S[, t] + + waning + rbinom(nlocs, size = S[, t-1] + V[, t-1], prob = births) 
    
    ##' Only do it all if vaccinated is greater than 0
    if (t %in% vstep) {
      n_vacc <- sim.vacc(locid = data$villcode, Sdogs = S[, t], Vdogs = V[, t],
                         vill_vacc = vacc[, t], p_revacc = p_revacc, p_allocate = p_vacc, 
                         type = vacc_sim, cov_est = cov_val, perc = perc_val)
    } else {
      n_vacc <- 0
    }
    
    ## balance vaccinated
    S[, t] <- S[, t] - n_vacc 
    V[, t] <- V[, t] + n_vacc
    
    # 3c. Exposed to infectious --------------------------------------------------------------------
    ## exposed to infectious
    E_coords$infectious <- rbinom(nrow(E_coords), size = 1, prob = sigma)
    
    ## add any that became infectious to I_coords_now
    I_coords_now <- E_coords[infectious == 1]
    
    ## remove any that became infectious from E_coords
    E_coords <- E_coords[infectious == 0] 
    
    ## Only count exposed -> infectious from within to subtract out
    if(!is.null(I_coords_now)) {
      if(nrow(I_coords_now) > 0) {
        I_count <- I_coords_now[, list(count = length(ID)), by = cell_id]
        I_dist[row_id[match(I_count$cell_id, cell_id)], t] <- I_count$count
      }
    } 
    
    # 3d. Add in incursions --------------------------------------------------------------------
    ## incursions (To do? : could add in option to weight by distance from edge)
    incursions <- rpois(1, inc_rate) # number of incursions in week
    
    if (incursions > 0){
      incs <- sim.incursions(incursions = incursions, counter = counter, row_ids = rows_pop, 
                             cell_ids = cell_id, tstep = t, x_coord = data$x_coord,
                             y_coord = data$y_coord)
      counter <- incs[["counter"]]
      I_locs <- incs[["I_locs"]]
      I_all[I_locs, t] <- 1
      I_coords_now <- rbindlist(list(incs[["I_coords_out"]], I_coords_now), fill = TRUE, use.names = TRUE)
    }
    
    I_all[, t] <- I_all[, t] + I_dist[, t]
    E[, t] <- E[, t-1] - I_dist[, t] ## - infectious from within district
    N[, t] <- S[, t] + V[, t] + E[, t] + I_dist[, t] ## everyone right now!
    
    # 3e. Simulate individual based contacts + movement --------------------------------------------
    E_new <- rep(0, nlocs)
    
    if(!is.null(I_coords_now)) {
      if(nrow(I_coords_now) > 0) {
        I_coords_now$secondaries <- rnbinom(nrow(I_coords_now), mu = R_0, size = k)
        ## Exposed in this time step
        exps <- sim.bites(secondaries = I_coords_now$secondaries, x_coord = I_coords_now$x_coord, 
                          y_coord = I_coords_now$y_coord, ids = I_coords_now$ID,
                          grid = grid, cells_pop = cells_pop, counter = counter, 
                          dispersalShape = dispersalShape, dispersalScale = dispersalScale, 
                          x_topl = x_topl, y_topl = y_topl, res_m = res_m, tstep = t)
        E_coords_now <- exps[["E_coords_now"]]
        counter <- exps[["counter"]]
        
        # 3f. Aggregate probability that contact is with susceptible -----------------------------------
        if(nrow(E_coords_now > 0)) {
          E_coords_now <- sim.trans(E_coords_now, cell_id = cell_id, row_id = row_id, S = S[, t],
                                    N = N[, t])
          if (nrow(E_coords_now) > 0) {
            E_count <- E_coords_now[, list(count = length(ID)), by = cell_id]
            E_new[row_id[match(E_count$cell_id, cell_id)]] <- E_count$count ## add in newly exposed
            ## Using rbindlist to speed up 
            E_coords <- rbindlist(list(E_coords, E_coords_now), fill = TRUE, use.names = TRUE)
          }
        }
      }  
    }
    
    ## Finish balancing everyone out
    E[, t] <- E[, t] + E_new ## + newly exposed
    S[, t] <- S[, t] - E_new ## - newly exposed
    

  if(return_coords == TRUE) {
    if(!is.null(I_coords_now)) {
      if(nrow(I_coords_now) > 0) {
          ## alternatively = fill in list of I_coords dataframe 
          I_coords <- rbindlist(list(I_coords, I_coords_now), fill = TRUE, use.names = TRUE) ## to build trees!
        }
      }
    }
  }
  
  if (return_coords == TRUE) {
    return(list(N = N, S = S, E = E, I_all = I_all, I_dist = I_dist, I_coords = I_coords))
  } else {
    return(list(N = N, S = S, E = E, I_all = I_all, I_dist = I_dist))
  }
}
