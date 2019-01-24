## HELPER FUNCTIONS FOR IBM
## MALAVIKA RAJEEV

# Transmission: first time step ---------------------------------------------------------------
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



# Transmission: subsequent timesteps ----------------------------------------------------------





# Allocating vacc -----------------------------------------------------------------------------


now.vill <- as.data.table(list(villcode = data$villcode, 
                               Sdogs = S[, t-1], Vdogs = V[, t-1], 
                               Ndogs = N[, t-1]))
now.vill <- now.vill[, .(Sdogs = sum(Sdogs, na.rm = TRUE), 
                         Vdogs = sum(Vdogs, na.rm = TRUE),
                         Ndogs = sum(Ndogs, na.rm = TRUE)), 
                     by = villcode]
now.vill$vacc <- vacc[, t + 1][match(now.vill$villcode, vacc$villcode)]
now.vill$revacc <- rbinom(nrow(now.vill), size = now.vill$Vdogs, prob = p_revacc)
## dogs currently vaccinated that were revaccinated
now.vill$vacc_new <- ifelse(now.vill$vacc - now.vill$revacc <= 0, now.vill$vacc, 
                            ifelse((now.vill$vacc - now.vill$revacc) > now.vill$Sdogs, 
                                   now.vill$Sdogs, 
                                   now.vill$vacc - now.vill$revacc))
## dogs newly vaccinated accounting for revacc(max)
## if revacc is bigger than # vaccinated total, likely indicates an additive, smaller campaign
## this is based on the way we grouped campaigns 
## TO DO: check to see what not assuming this does to vacc cov
## otherwise just maximum possible for that vill

## Keep track of village cov
if (t == 2){
  rownames(vill_vacc) <- now.vill$villcode
}

vill_vacc[, t] <- (now.vill$Vdogs + now.vill$vacc_new)/now.vill$Ndogs

## Get cell level vaccinated
## Get vacc prob for that grid proportional to total pop in that vill
## basically means that places with more dogs are more likely to be vaccinated
## Could also do it by the available # of dogs (sus_grid/sus_vill)
vacc_now <- data[now.vill, on = "villcode"]
## make sure order is preserved!  
vacc_now <- vacc_now[order(match(vacc_now$cell_id, data$cell_id)), ]
vacc_now$sus <- S[, t-1]
vacc_now[, vacc_prob := ifelse(Sdogs == 0, 1e-6, sus/Sdogs)]
vacc_now[, n_vacc := as.double(rmultinom(1, size = first(vacc_new), 
                                         prob = vacc_prob)), 
         by = villcode]

## Set so that newly vaccinated cannot exceed the number sus in that grid cell
## Way to do 1e-6 multinomial cheat and still be okay 
## Keeping track of vacc individuals that we lose
vacc_now$leftovers <- ifelse(vacc_now$n_vacc > vacc_now$sus, 
                             vacc_now$n_vacc - vacc_now$sus, 0)
vacc_now$n_vacc <- ifelse(vacc_now$n_vacc > vacc_now$sus, vacc_now$sus,
                          vacc_now$n_vacc)

## need to change unnacounted so that it's about the vill level sum instead!
unaccounted[t] <- sum(vacc_now$leftovers, na.rm = TRUE)
