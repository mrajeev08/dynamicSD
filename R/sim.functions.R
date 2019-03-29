##' 1. Simulate vaccination ------------------------------------------------------------------------
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
#' 
#' @return Vector of either same length as number of grid cells with number of newly vaccinated dogs 
#'   in each grid cell or 0L returned when no vaccination occured

sim.vacc <- function(locid = data$villcode, Sdogs, Ndogs, Vdogs, p_revacc, p_allocate = 1, 
                     vill_vacc, ...) {
  
  ##' TO DO:
  ##' checks on 'lost vaccinations' at village level
  ##' check to see what assumptions about revaccination and additive campaigns do to cov estimates
  ##' compare allocation of vacc to see what diff it makes

  vacc_grid <- as.data.table(list(locid = locid, Sdogs = Sdogs, Vdogs = Vdogs, Ndogs = Ndogs, 
                                  vacc = vill_vacc))
  vacc_vill <- vacc_grid[, .(Sdogs = sum(Sdogs, na.rm = TRUE), 
                           Vdogs = sum(Vdogs, na.rm = TRUE),
                           Ndogs = sum(Ndogs, na.rm = TRUE), 
                           vacc = first(vacc)), 
                       by = locid]
  
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
    vacc_sample <- as.data.table(list(order_id = rep(vacc_grid$order_id, times = vacc_grid$sus)))
    vacc_sample <- vacc_sample[vacc_grid, on = "order_id"]
    vaccs <- vacc_sample[, .(order_id = sample(order_id, size = first(vacc_new), 
                                            replace = FALSE, prob = p_allocate)),
                         by = "locid"]
    vaccs <- vaccs[,.(n_vacc = .N), by = "order_id"]
    vacc_grid <- vacc_grid[vaccs, on = "order_id"]
    vacc_grid$n_vacc[is.na(vacc_grid$n_vacc)] <- 0
  }

  
  ##' Set so that newly vaccinated cannot exceed the number sus in that grid cell
  ##' If we set p_allocate to the prob from the census, then we might lose lots of potential 
  ##' vaccinations this way
  vacc_grid[, n_vacc := as.double(rmultinom(1, size = first(vacc_new), 
                                            prob = vacc_prob)), by = locid]
  vacc_grid$n_vacc <- ifelse(vacc_grid$n_vacc > vacc_grid$sus, vacc_grid$sus,
                             vacc_grid$n_vacc)
  ##' To do = add as test!
  ##' Should end up in same order so that (try length(vacc_grid$locid == locid)[FALSE], should = 0)
  
  n_vacc <- vacc_grid$n_vacc
  return(n_vacc)
}
