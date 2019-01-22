## Function for running batches of sims
sim.IBM <- function(nsims = 10, grid = SD_raster, data = grid_data, vacc = vacc_mat,
                    I_seed = 2, start_vacc = 0.2,  
                    R_0 = 1.1, k = 0.2, 
                    inc_rate = incs,
                    # weekly number of incursions and scaling factor 
                    sigma = get.prob(rate = 7/22.3, step = 1), # weekly rate to prob exp -> inf
                    births = get.prob(rate = grid_data$births, step = 52),
                    cell_id = grid_data$cell_id,
                    nlocs = nrow(grid_data),
                    # annual birth rate to prob
                    mu = get.prob(rate = 0.44, step = 52), # annual death rate to prob 
                    ntimes = tmax, # maximum time step,
                    nu = get.prob(rate = 0.33, step = 52), p_revacc = 0.5, 
                    # annual waning rate to prob + probability of revaccination
                    dispersalShape = 0.3484, dispersalScale = 41.28/100, # from Rebecca
                    return_coords = TRUE) {
  ## Summary matrices
  row_id <- 1:nlocs
  I_mat_ts <- S_mat_ts <- N_mat_ts <- matrix(NA, nrow = ntimes, ncol = nsim)
  I_mat_loc <- S_mat_loc <- N_mat_loc <- matrix(NA, nrow = ntimes, ncol = nsim)
  
  foreach(i = 1:nsims)
  ) %dopar% { 
    # 1. Set-up and tstep 1 -----------------------------------------------------------------------
    
    # state matrices
    S <- V <- N <- matrix(0, nrow = nlocs, ncol = ntimes)
    I_dist <- I_all <- E <- matrix(0, nrow = nlocs, ncol = ntimes)
    
    ## Starting pop + sus 
    N[, 1] <- data$start_pop
    V[, 1] <- rbinom(n = nlocs, size = data$start_pop, prob = start_vacc)
    S[, 1] <- N[, 1] - V[, 1]
    
    # 2. Simulating infection in first timestep ----------------------------------------------------
    # Seeding cases (as if from outside!)
    # pick random location to start I seeds
    I_dist[ ,1] <- 0
    I_locs <- sample((row_id)[which(!is.na(data$start_pop))], I_seed)
    I_dist[I_locs, 1] <- 1
    
    ## Create infected coordinate data frame to grow
    I_coords <- data.table(ID = 1:I_seed, tstep = rep(1, I_seed), x_coord = data$x_coord[I_locs], 
                           y_coord = data$y_coord[I_locs], progen_ID = rep(0, I_seed),
                           path_ID = rep(NA, I_seed), cell_id = cell_id[I_locs], 
                           sus = NA, trans = NA,
                           infectious = rep(1, I_seed), secondaries = NA)
    
    ## Create exposed coordinate data frame to grow
    E_coords <- data.table(ID = NA, tstep = NA, x_coord = NA, y_coord = NA, 
                           progen_ID = NA, path_ID = NA, cell_id = NA, sus = NA, trans = NA, 
                           infectious = NA, secondaries = NA)
    
    counter <- max(I_coords$ID) ## for keeping track of ids
    
    if(nrow(I_coords) > 0){
      I_coords$secondaries <- rnbinom(nrow(I_coords), mu = R_0, size = k)
      for (i in 1:nrow(I_coords)){
        if (I_coords$secondaries[i] > 0) {
          for (j in 1:I_coords$secondaries[i]){ 
            counter <- counter + 1
            if (j == 1) { # need progenitor coords for 1st movement
              origin_x <- I_coords$x_coord[i]
              origin_y <- I_coords$y_coord[i]
            } else {
              origin_x <- last_coords[1]
              origin_y <- last_coords[2]
            }
            within <- 0
            while (within == 0) {
              distance <- rweibull(1, shape = dispersalShape, 
                                   scale = dispersalScale) # Distance to move in km
              angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at
              
              # Convert to a vector and move
              x_new <- (sin(angle) * distance * 1000) + origin_x # convert to m
              y_new <- (cos(angle) * distance * 1000) + origin_y
              
              ## This means will only go somewhere within a vill
              check <- extract(grid, cbind(x_new, y_new))
              if(check %in% cell_id){
                within <- ifelse(N[row_id[which(cell_id == check)], 1] > 0, 1, 0)
              }
              ## This code means only in populated places and doggies don't leave the district
              ## or not restricting to populated places:
              ## within <- ifelse(check %in% data$cell_id == TRUE, 1, 0)
            }
            
            E_coords_new <- list(ID = counter, tstep = 1, 
                                 x_coord = x_new, y_coord = y_new,  
                                 progen_ID = I_coords$ID[i], path_ID = j, 
                                 cell_id = check, sus = NA, trans = 0, 
                                 infectious = NA, secondaries = NA)
            last_coords <- c(x_new, y_new)
            E_coords <- rbindlist(list(E_coords, E_coords_new), fill = TRUE, use.names = TRUE)
          }
        }
      }
    }
    
    E_coords <- E_coords[!is.na(ID)]
    
    if(nrow(E_coords) > 0) {
      ## probability that contact will be with a susceptible = St/Nt
      E_coords$sus <- S[row_id[match(E_coords$cell_id, cell_id)], 1]
      E_coords$N <- N[row_id[match(E_coords$cell_id, cell_id)], 1]
      ids <-  unique(E_coords$cell_id)
      for (i in 1:length(ids)) {
        E_trans <- E_coords[cell_id == ids[i]]
        done <- 0
        max <- E_trans$sus[1]
        for(j in 1:nrow(E_trans)){
          if (done < max){
            E_trans$trans[j] <- rbinom(1, size = 1, prob = ifelse(E_trans$N[j] == 0, 0, 
                                                                  E_trans$sus[j]/E_trans$N[j]))
            E_coords$trans[match(E_trans$ID[j], E_coords$ID)] <- E_trans$trans[j]
            done <- done + E_trans$trans[j]
          }
        }
      }
      E_coords$sus <- E_coords$sus/E_coords$N
      E_coords[, N := NULL]
      E_coords <- E_coords[!is.na(trans) & trans == 1]
    }
    
    if (nrow(E_coords) > 0) {
      E_count <- E_coords[, list(count = length(ID)), by = cell_id]
      E[row_id[match(E_count$cell_id, cell_id)], 1] <- E_count$count ## add in newly exposed
    }
    
    ## Finish subtracting everyone out
    S[, 1] <- S[, 1] - E[ ,1] ## subtract out newly exposed guys from S
    unaccounted <- rep(0, ntimes) ## keeping track of any vacc individuals that went unaccounted for
    
    # 3. Simulate for rest of time steps -----------------------------------------------------------
    for (t in 2:ntimes) {
      # t = 2
      # print(paste(t, "/", tmax, "weeks"))
      
      # 3a. Simulate vaccination -------------------------------------------------------------------
      ## TO DO: ## make sure leftover vaccinated goes to adjacent vills
      now.vill <- as.data.table(list(villcode = data$villcode, 
                                     sus = S[, t-1], vdogs = V[, t -1]))
      now.vill <- now.vill[, .(sus = sum(sus, na.rm = TRUE), vdogs = sum(vdogs, na.rm = TRUE)), 
                           by = villcode]
      now.vill$vacc <- vacc[, t + 1][match(now.vill$villcode, vacc$villcode)]
      now.vill$revacc <- rbinom(nrow(now.vill), size = now.vill$vdogs, prob = p_revacc)
      now.vill$vacc_new <- ifelse(now.vill$vacc - now.vill$revacc <= 0, now.vill$vacc, 
                                  ifelse((now.vill$vacc - now.vill$revacc) > now.vill$sus, 
                                         now.vill$sus, 
                                         now.vill$vacc - now.vill$revacc)) 
      ## if revacc is bigger than # vaccinated total, likely indicates an additive campaign
      ## otherwise just maximum possible for that vill
      data$vacc_new <- now.vill$vacc_new[match(data$villcode, now.vill$villcode)]
      data$sus <- S[, t-1]
      
      ## trying to use data table instead
      vacc_now <- data[, vacc_prob := ifelse(vacc_new == 0, 0, sus/vacc_new), by = villcode]
      vacc_now[, n_vacc := as.double(ifelse(vacc_new == 0, 0, rmultinom(1, size = first(vacc_new), 
                                                                        prob = vacc_prob))), by = villcode]
      vacc_now$n_vacc <- ifelse(vacc_now$n_vacc - vacc_now$sus > 0, vacc_now$sus, vacc_now$n_vacc)
      vacc_now$leftovers <- ifelse(vacc_now$n_vacc - vacc_now$sus > 0, 
                                   vacc_now$n_vacc - vacc_now$sus, 0)
      unaccounted[t] <- sum(vacc_now$leftovers, na.rm = TRUE)
      ## need to change unnacounted so that it's about the vill level sum instead!
      
      ## ends up in same order (try length(vacc_now$cell_ID == data$cell_ID)[FALSE], should = 0)
      
      # 3b. Sequential transitions for S + V -----------------------------------------------
      ## S transitions, sequenctially for competing probs (-vacc, - deaths)
      S[, t] <- S[, t-1] - vacc_now$n_vacc ## - vaccinated first
      S[, t] <- S[, t] - rbinom(nlocs, size = S[, t], prob = mu) ## then die
      
      ## V transitions, sequetially for competing probs (-waning, - deaths)
      waning <- rbinom(nlocs, size = V[, t-1], prob = nu)
      V[, t] <- V[, t-1] - waning
      V[, t] <- V[, t] - rbinom(nlocs, size = V[, t], prob = mu)
      
      ## Additive probs
      V[, t] <- V[, t] + vacc_now$n_vacc
      S[, t] <- S[, t] + rbinom(nlocs, size = S[, t-1] + V[, t-1], prob = births) + waning
      
      # 3c. Exposed to infectious --------------------------------------------------------------------
      ## exposed to infectious
      E_coords$infectious <- rbinom(nrow(E_coords), size = 1, prob = sigma)
      
      I_coords_in <- E_coords[!is.na(infectious) & infectious == 1]
      ## add any that became infectious to I_coords_now
      
      ## Only count exposed -> infectious from within to subtract out
      if(nrow(I_coords_in) > 0) {
        I_count <- I_coords_in[, list(count = length(ID)), by = cell_id]
        I_dist[row_id[match(I_count$cell_id, cell_id)], t] <- I_count$count
      } 
      
      ## remove any that became infectious from E_coords
      E_coords <- E_coords[!is.na(infectious) & infectious == 0] 
      
      # 3d. Add in incursions --------------------------------------------------------------------
      ## incursions (To do? : could add in option to weight by distance from edge)
      incursions <- rpois(1, inc_rate) # number of incursions in week
      if (incursions > 0){
        start_counter <- counter + 1 
        counter <- counter + incursions
        I_all[, t] <- 0
        
        ## only happen in populated places
        I_locs <- sample((row_id)[which(!is.na(data$start_pop))], incursions) 
        I_all[I_locs, t] <- 1 
        
        I_coords_out <- data.table(ID = start_counter:counter, tstep = t, 
                                   x_coord = data$x_coord[I_locs],
                                   y_coord = data$y_coord[I_locs], progen_ID = 0, path_ID = NA, 
                                   cell_id = data$cell_id[I_locs], sus = NA,
                                   trans = NA, infectious = 1, 
                                   secondaries = NA)
        I_coords_now <- rbindlist(list(I_coords_in, I_coords_out), fill = TRUE, use.names = TRUE)
      } else {
        I_coords_now <- I_coords_in
      }
      
      I_all[, t] <-I_all[, t] + I_dist[, t]
      
      # 3e. Simulate individual based contacts + movement --------------------------------------------
      if(nrow(I_coords_now) > 0){
        E_coords_now <- E_coords[nrow(E_coords)+1] ## to get empty df to bind to!
        I_coords_now$secondaries <- rnbinom(nrow(I_coords_now), mu = R_0, size = k)
        
        for (i in 1:nrow(I_coords_now)){
          if (I_coords_now$secondaries[i] > 0) {
            for (j in 1:I_coords_now$secondaries[i]){ 
              counter <- counter + 1
              if (j == 1) { # need progenitor coords for 1st movement
                origin_x <- I_coords_now$x_coord[i]
                origin_y <- I_coords_now$y_coord[i]
              } else {
                origin_x <- last_coords[1]
                origin_y <- last_coords[2]
              }
              within <- 0
              while (within == 0) {
                distance <- rweibull(1, shape = dispersalShape, 
                                     scale = dispersalScale) # Distance to move in km
                angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at
                
                # Convert to a vector and move
                x_new <- (sin(angle) * distance * 1000) + origin_x # convert to m
                y_new <- (cos(angle) * distance * 1000) + origin_y
                
                ## This means will only go somewhere within a vill
                check <- extract(grid, cbind(x_new, y_new))
                
                if(check %in% cell_id){
                  within <- ifelse(N[row_id[which(cell_id == check)], t-1] > 0, 1, 0)
                }
                ## This code means only in populated places and doggies don't leave the district
                ## or not restricting to populated places:
                ## within <- ifelse(check %in% data$cell_id == TRUE, 1, 0)
              }
              
              E_coords_new <- list(ID = counter, tstep = t, 
                                   x_coord = x_new, y_coord = y_new,  
                                   progen_ID = I_coords_now$ID[i], path_ID = j, 
                                   cell_id = check, sus = NA, trans = NA, 
                                   infectious = 0, secondaries = NA)
              last_coords <- c(x_new, y_new)
              E_coords_now <- rbindlist(list(E_coords_now, E_coords_new), fill = TRUE, use.names = TRUE)
            }
          }
        }
        E_coords_now <- E_coords_now[-1, ] ## take out NA row generated to bind to df
      } else {
        E_coords_now <- I_coords_now ## gets you empty df
      }
      
      # 3f. Aggregate probability that contact is with susceptible -----------------------------------
      E[, t] <- E[, t-1] - I_dist[, t] ## - infectious from within district + newly exposed
      N[, t] <- S[, t] + V[, t] + E[, t] + I_dist[, t]
      
      if(nrow(E_coords_now) > 0) {
        ## probability that contact will be with a susceptible = St/Nt-1
        E_coords_now$sus <- S[row_id[match(E_coords_now$cell_id, cell_id)], t]
        E_coords_now$N <- N[row_id[match(E_coords_now$cell_id, cell_id)], t]
        ids <-  unique(E_coords_now$cell_id)
        setkey(E_coords_now, cell_id)
        for (i in 1:length(ids)) {
          E_trans <- E_coords_now[J(ids[i])]
          done <- 0
          max <- E_trans$sus[1]
          for(j in 1:nrow(E_trans)){
            if (done < max){
              E_trans$trans[j] <- rbinom(1, size = 1, prob = ifelse(E_trans$N[j] == 0, 0, 
                                                                    E_trans$sus[j]/E_trans$N[j]))
              E_coords_now$trans[match(E_trans$ID[j], E_coords_now$ID)] <- E_trans$trans[j]
              done <- done + E_trans$trans[j]
            }
          }
        }
        
        E_coords_now$sus <- E_coords_now$sus/E_coords_now$N
        E_coords_now[, N := NULL]
        E_coords_now <- E_coords_now[!is.na(trans) & trans == 1]
        
      } else {
        E_coords_now <- E_coords_now # to get empty df to bind!
      }
      
      E_new <- rep(0, nlocs)
      
      if (nrow(E_coords_now) > 0) {
        E_count <- E_coords_now[, list(count = length(ID)), by = cell_id]
        E_new[row_id[match(E_count$cell_id, cell_id)]] <- E_count$count ## add in newly exposed
      }
      
      ## Finish subtracting everyone out
      E[, t] <- E[, t] + E_new ## - infectious from within district + newly exposed
      S[, t] <- S[, t] - E_new ## subtract out newly exposed guys in the week
      
      ## change this to avoid growing dataframes!
      E_coords <- rbindlist(list(E_coords, E_coords_now), fill = TRUE, use.names = TRUE)
      I_coords <- rbindlist(list(I_coords, I_coords_now), fill = TRUE, use.names = TRUE) ## to build trees!
    }
    if (return_coords == TRUE) {
      return(list(N, S, E, I_all, I_dist, I_coords))
    }
    if (return_coords == FALSE) {
      return(list(N, S, E, I_all, I_dist))
    }
  }
}