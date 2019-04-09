##' Simulate vaccination 
##' ------------------------------------------------------------------------------------------------
#' Simulate vaccination campaigns spatially 
#'
#' \code{sim.vacc} returns the dogs newly vaccinated in each grid cell. 
#'
#' This function takes the input of dogs vaccinated (weekly or monthly, at whatever scale you aggregate)
#'   in each location(i.e. village for Serengeti) and allocates them to grid cells within a village 
#'   accounting for revaccination using multinomial or sampling by probabilities. 
#'   It is critical that all vectors are in same order.
#'
#' @param locid the location id to match between grid cells and data in 'vill_vacc' 
#' @param Sdogs number of susceptible dogs in each grid cell at time t
#' @param Vdogs number of vaccinated dogs in each grid cell at time t
#' @param Ndogs number of total dogs in each grid cell at time t (includes exposed)
#' @param p_revacc probability of revaccination of currently vaccinated dogs
#' @param p_allocate vector of probabilities of vaccination for each grid cell (currently this can either
#'   be by the proportion of the population in each grid cell at the village level or 
#'   the reported vaccination coverage in that grid cell estimated from the HH census. 
#'   @seealso [get.demdata()]). If set to 1 then vaccinations are allocated by the proportion of 
#'   susceptible hosts available to be vaccinated at the village level that are in the grid cell. 
#' @param vill_vacc vector of # vaccinated by village at time t for each grid cell
#' @param type vacc data being passed through, either simulated at cell level ("sim_cell"),
#'   simulated at village level ("sim_vill"), or empirical # vaccinated from SD ("empirical")
#' @return Vector of either same length as number of grid cells with number of newly vaccinated dogs 
#'   in each grid cell or 0L returned when no vaccination occured
#'   
#' @section Dependencies:
#'     Packages: data.table, 

sim.vacc <- function(locid = data$villcode, Sdogs, Ndogs, Vdogs, p_revacc, p_allocate = 1, 
                     vill_vacc, type = "empirical", cov_est = 1, perc = 1,...) {
  
  ##' TO DO:
  ##' checks on 'lost vaccinations' at village level
  ##' check to see what assumptions about revaccination and additive campaigns do to cov estimates
  ##' compare allocation of vacc to see what diff it make
  ##' 
  if(type == "sim_cell") {
    vill_vacc <- rbinom(1, size = 1, prob = perc*vill_vacc) 
    ## what's the prob that it will get vaccinated at given cov level
    vacc <- rbinom(length(Ndogs), size = Ndogs, prob = vill_vacc*cov_est) 
    ## if vaccinated in grid cell, what the cov est abt
    
    revacc <- rbinom(length(Vdogs), size = Vdogs, prob = p_revacc)
    n_vacc <- ifelse(vacc - revacc <= 0, 0, 
                     ## if higher est of cov then vacc_vill$vacc, if mid est, then vacc_vill$vacc*p_revacc
                     ifelse((vacc - revacc) > Sdogs, Sdogs, 
                            vacc - revacc))
    
  }
  
  else {
    vacc_grid <- as.data.table(list(locid = locid, Sdogs = Sdogs, Vdogs = Vdogs, Ndogs = Ndogs, 
                                    vacc = vill_vacc))
    vacc_vill <- vacc_grid[, .(Sdogs = sum(Sdogs, na.rm = TRUE), 
                             Vdogs = sum(Vdogs, na.rm = TRUE),
                             Ndogs = sum(Ndogs, na.rm = TRUE), 
                             vacc = first(vacc)), 
                         by = locid]
    
    ## if working with simulated campaigns (not specifying # of dogs vaccinated)
    if(type == "sim_vill") {
      vacc_vill$vacc <- rbinom(1, size = 1, prob = perc*vacc_vill) 
      vacc_vill[, vacc := rbinom(1, size = Ndogs, prob = vacc*cov_est)]
    } 
    
    ##' dogs currently vaccinated that were revaccinated
    vacc_vill$revacc <- rbinom(nrow(vacc_vill), size = vacc_vill$Vdogs, prob = p_revacc)
      
    ##' dogs newly vaccinated accounting for revacc(max)
    ##' if revacc is bigger than # vaccinated total, likely indicates an additive, smaller campaign?
    vacc_vill$vacc_new <- ifelse(vacc_vill$vacc - vacc_vill$revacc <= 0, 0, 
                                 ## if higher est of cov then vacc_vill$vacc, if mid est, then vacc_vill$vacc*p_revacc
                                 ifelse((vacc_vill$vacc - vacc_vill$revacc) > vacc_vill$Sdogs, 
                                        vacc_vill$Sdogs, 
                                        vacc_vill$vacc - vacc_vill$revacc))
    
    ##' Get cell level vaccinated (need to do it this way so that susceptibles don't get lost)
    vacc_grid <- as.data.table(list(locid = locid, sus = Sdogs, p_allocate = p_allocate,
                                    order_id = 1:length(Sdogs)))
    vacc_grid <- vacc_grid[vacc_vill, on = "locid"]
    
    ##' make sure original order of df is preserved!  
    vacc_grid <- vacc_grid[order(vacc_grid$order_id), ]   
    
    ##' Allocate by available susceptibles at the current time in the village and constrain to 1e-10
    ##'  so that multinomial will still work
    ##' basically means that places with more susceptible dogs are more likely to be vaccinated
    if(length(p_allocate) == 1) {
      vacc_grid[, vacc_prob := ifelse(Sdogs == 0 , 1e-10, sus/Sdogs)]
      vacc_grid[, n_vacc := as.double(rmultinom(1, size = first(vacc_new), prob = vacc_prob)), 
                by = locid]
    ##' Or allocate by other vaccination probability (this will take longer for sure)
    } else {
      vacc_grid$sus[is.na(vacc_grid$sus)] <- 0
      vacc_sample <- as.data.table(list(order_id = rep(vacc_grid$order_id, times = vacc_grid$sus)))
      vacc_sample <- vacc_sample[vacc_grid, on = "order_id"]
      vaccs <- vacc_sample[, .(order_id = sample(order_id, size = first(vacc_new), 
                                              replace = FALSE, prob = p_allocate)),
                           by = "locid"]
      vaccs <- vaccs[,.(n_vacc = .N), by = "order_id"]
      vacc_grid <- vaccs[vacc_grid, on = "order_id"]
      vacc_grid$n_vacc[is.na(vacc_grid$n_vacc)] <- 0
    }
    
    ##' Set so that newly vaccinated cannot exceed the number sus in that grid cell
    vacc_grid$n_vacc <- ifelse(vacc_grid$n_vacc > vacc_grid$sus, vacc_grid$sus,
                               vacc_grid$n_vacc)
    ##' To do = add as test!
    ##' Should end up in same order so that (try length(vacc_grid$locid == locid)[FALSE], should = 0)
    
    n_vacc <- vacc_grid$n_vacc
  }
  return(n_vacc)
}

##' Simulating incursions 
##' ------------------------------------------------------------------------------------------------
#' Simulate incursions spatially
#' 
#' \code{sim.incursions} simulates incursions draw from a avg number per timestep and assigns spatially.
#' 
#' This function 
#' @param Paramters
#' @return Returned
#' @section Dependencies:
#'     List dependencies here, i.e. packages and other functions

sim.incursions <- function(incursions, counter, row_ids, cell_ids, tstep, 
                           x_coord, y_coord) {
  start_counter <- counter + 1 
  counter <- counter + incursions
  I_all <- rep(0, length(cell_ids))
    
  ## only happen in populated places
  I_locs <- sample(row_ids, incursions) 
    
  I_coords_out <- data.table(ID = start_counter:counter, tstep = tstep, 
                               x_coord = x_coord[I_locs],
                               y_coord = y_coord[I_locs], progen_ID = 0, path_ID = NA, sus = NA,
                               trans = NA, cell_id = cell_ids[I_locs], infectious = 1, 
                               secondaries = NA)
  return(list(I_coords_out = I_coords_out, I_locs = I_locs, counter = counter))
}



##' Simulate biting
##'  -----------------------------------------------------------------------------------------------
#' Simulate biting and movement
#' 
#' \code{sim.bites} simulates bites and movement on landscape
#' 
#' Number of secondary cases is drawn from a negative binomial distribution with mean R0 and 
#' dispersion parameter k. Movement is sequential so that the dispersal kernal is movement 
#' between bite events (as per Rebecca). Cell ids are determined by subtracting 
#' from the top left corner of grid to get 1-based indexes (use bbox to figure this out).
#' 
#' @param secondaries number of secondary cases for each infectious at tstep
#' @param x_coord UTM Easting for each infectious case at tstep
#' @param y_coord UTM Northing for each infectious case at tstep
#' @param ids case IDs for each infectious case at tstep
#' @param counter the counter for the ID of each potential infectious individual
#' @param grid matrix of cell_ids from gridded raster of study area 
#' @param res_m resolution of raster in meters 
#' @param cells_pop populated cells in raster (i.e. where dogs are allowed to move, could also define
#'   this in other ways)
#' @param dispersalShape shape parameter of Weibull parameterized as per Rebecca.
#' @param dispersalScale scale parameter of Weibull 
#' @param x_topl top left x coordinate of grid, to calculate cell_id based on 1-based cell indexing
#' @param y_topl top left y coordinate of grid
#' @param tstep current timestep 
#' 
#' @return Returns data.table of exposures in current timestep with following columns: ID, tstep, 
#'   x_coord, y_coord, progen_ID, path_ID, cell_id, sus, trans, infectious, secondaries, 
#'   
#' @section Dependencies:
#'     Packages: None

sim.bites <- function(secondaries = I_coords_now$secondaries, x_coord = I_coords_now$x_coord, 
                      y_coord = I_coords_now$y_coord, ids = I_coords_now$ID, counter, grid, res_m, 
                      cells_pop, dispersalShape, dispersalScale, x_topl, y_topl, tstep = t) {
  ##' TO DO:
  ##' Vectorize movements between progenitors (i.e. the number of movements = # of progenitors)
  ##' Only really need to do that if dealing with lots of infected in one timestep
  ##' 
  ##' Need to think about what makes sense in terms of restricting movement
  ## if not restricting to populated places then cells should be anywhere in the district
  ##' if not restricting to district then doggies can leave if they want so then either
  ##' they bite in a populated place within the district or they go outside the district 
  ##' which means coords are in greater than x_min/y_max zone...
  ##' ...need to think about what's plausible and whether it would make a difference... 
  
  store_coords_all <- vector("list", length(secondaries)) 
  
  for (i in 1:length(secondaries)) {
    if (secondaries[i] > 0) {
      
      store_coords <- vector("list", secondaries[i]) 
      
      for (j in 1:secondaries[i]){ 
        counter <- counter + 1
        if (j == 1) { # need progenitor coords for 1st movement
          origin_x <- x_coord[i]
          origin_y <- y_coord[i]
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
          
          ## This means can go anywhere including outside of district
          col <- ceiling((x_new - x_topl)/res_m)
          row <- ceiling(-(y_new - y_topl)/res_m)
          cell <- ifelse(row %in% 1:nrow(grid) & col %in% 1:ncol(grid), grid[row, col], 0)
          ## This bit means only in populated places and doggies don't leave the district
          
          within <- ifelse(cell %in% cells_pop, 1, 0)
        }
  
        store_coords[[j]] <- list(ID = counter, tstep = tstep, x_coord = x_new, y_coord = y_new,  
                                  progen_ID = ids[i], path_ID = j, cell_id = cell, sus = NA, 
                                  trans = NA_integer_, infectious = 0L, secondaries = NA_integer_)
        last_coords <- c(x_new, y_new)
      }
      store_coords_all[[i]] <- store_coords
    }
  }
  E_coords_now <- rbindlist(lapply(store_coords_all, rbindlist))
  return(list(E_coords_now = E_coords_now, counter = counter))
}

##' Simulate transmission
##'  -----------------------------------------------------------------------------------------------
#' Simulating transmission 
#' 
#' \code{sim.trans} simulates transmission events (i.e. whether a contact (bite) is with a 
#' susceptible individual)
#' 
#' Details
#' 
#' @param Parameters
#' @return Returned
#' @section Dependencies:
#'     List dependencies here, i.e. packages and other functions

sim.trans <- function(E_coords_now, cell_id, row_id, S, N) {

  ## probability that contact will be with a susceptible = St/Nt
  E_coords_now$sus <- S[row_id[match(E_coords_now$cell_id, cell_id)]]
  E_coords_now$N <- N[row_id[match(E_coords_now$cell_id, cell_id)]]

  ids <-  unique(E_coords_now$cell_id)
  setkey(E_coords_now, cell_id)

  for (i in 1:length(ids)) {
    E_trans <- E_coords_now[J(ids[i])]
    sus <- E_trans$sus[1]
    
    for(j in 1:nrow(E_trans)){
      if (sus > 0){
        trans <- rbinom(1, size = 1, prob = ifelse(E_trans$N[j] == 0, 0,
                                                   sus/E_trans$N[j]))
        sus <- sus - trans
        E_coords_now$trans[match(E_trans$ID[j], E_coords_now$ID)] <- trans
      }
    }
  }
  
  E_coords_now$sus <- E_coords_now$sus/E_coords_now$N
  E_coords_now[, N := NULL]
  E_coords_now <- E_coords_now[!is.na(trans) & trans == 1]
  return(E_coords_now)
}


