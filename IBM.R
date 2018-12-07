### Functions for running most pared down IBM possible
## Malavika Rajeev
## December 2018

## Data
# shapefile
# dog census
# human census

## Step 1: get birth rate at village level
get.births <- function(grid, vill_vacc, pop02, pop12, census_pop)


## Step 2: Get vacc matrix at village level


## Step 3: Utility functions for IBM
## Get movement of infected individuals


## Simulate multinomial probabilities for each transition
sim.multinom <- function(sizes = c(10, 15, 20, 30), 
                         probs = c(0.01, 0.05)) {
  if((1 - sum(probs)) != 0) {
    prob_vec <- c(probs, 1 - sum(probs))
  } else {
    prob_vec <- probs
  }
  
  trans <- matrix(NA, length(prob_vec), length(sizes))
  
  for (i in 1:length(sizes)){
    trans[, i] <- rmultinom (1, size = sizes[i], prob_vec)
  }
  return(trans) # gets rid of other probabilties
}

sim.IBM <- function(shapefile, res_meters = 1000, N_0, vacc_mat, adj_mat
                    I.seed = 2, E.seed = 2, 
                    R_0 = 1.1, k = 0.2, 
                    sigma = 7/22.3,  
                    mu = (1 + 0.44)^(1/52)-1, br = births, 
                    tmax = 300, 
                    waning = (1 + 0.33)^(1/52)-1,
                    dispersalShape, dispersalScale) {
  
  ## check to see if it exists, if it doesn't make it!
  grid <- raster(shapefile)
  res(grid) <- res_meters
  grid <- rasterize (shapefile, grid)
  grid_data <- as.data.frame(grid)
  
  ## merge with pop + births and allocate multinomially across villages
  
  ## add in start_vacc
  
  ## vacc prob 
  nlocs <- nrow(grid_data)
  ntimes <- tmax + 1
  
  # state matrices
  S <- E <- I <- V <- N <- matrix(NA, nrow = nlocs, ncol = ntimes)
  
  # Starting indices 
  S[ ,1] <- start_sus
  N[ ,1] <- start_pop
  V[ ,1] <- start_pop - start_sus
  
  # Seeding cases
  # pick random location to start I/E seed 
  E[ ,1] <- 0
  E[sample(nlocs, Eseed), 1] <- 1
  
  I[ ,1] <- 0
  I[sample(nlocs, Iseed), ] <- 1
  S_probs <- c(births = br, )
  
  for (t in 2:ntimes){
    grid_data %>% group_by(vill) %>% mutate(pop = sum(N), prob_pop = pop/sum(N))
    left_join(grid_data, vacc_mat(indexed))
    mutate(prob_vacc = vacc_vill/pop, vacc = prob_vacc*prob_pop)
    
    S_trans <- sim.multinom(sizes = S[, t-1], prob = c(vacc, mu))
    E_trans <- sim.multinom(sizes = E[, t-1], prob = c(sigma))
    V_trans <- sim.multinom(sizes = V[, t-1], prob = c(waning, mu))
    E_trans[1, ] <- infection
    
    ## incursions
    
    ## draw movement for progenitors
    
    ## draw movement for secondary 
    for (i in 1:nlocs){
      exposures[i] = sum(rnbinom(n = infection[i] + incursions[i], mu = R0, size = k))
      for (j in 1:exposures[i]){
        # draw movement
        
        # assign to location 
        
        # keep track of ones who go outside
        # 
      }
    }
    
    transmission <- rbinom(nlocs, exposed, sus_prob)
    
    ## vaccination
    
    I[, t] <- infection
    E[, t] <- E[, t-1] - infection + transmission
    V[, t] <- V[, t-1] + vacc - V_trans[1, ] - V_trans[2, ]
    S[, t] <- S[, t-1] + S_trans[1, ] - S_trans[2, ] + V_trans[1, ] - vacc - transmission
  }
  ## return matrices
}             


GetMvtVector <- function(originX, originY, dispersalShape, dispersalScale, minX, maxX, minY, maxY){

    distance <- rweibull(1, shape=dispersalShape, scale=dispersalScale) # Distance to move in km
    angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at
    
    # Convert to a vector and move
    x.new <- (sin(angle) * distance * 1000) + originX # convert to m
    y.new <- (cos(angle) * distance * 1000) + originY
    
    # If within area keep, otherwise draw new vectors
    if(x.new>minX & x.new<maxX & y.new>minY & y.new<maxY){
      xy <- c(x.new, y.new); names(xy) <- c("x","y")
    }
  }
  return(xy)
}

## Testing
nrow(SD_coords)
Sys.time()
E_coords <- data.frame(ID = NA, tstep = NA, x_coord = NA, y_coord = NA, progen_ID = NA,
                       path_ID = NA)
I_coords <- data.frame(ID = NA, tstep = NA, x_coord = NA, y_coord = NA, progen_ID = NA)

SD_coords <- coordinates(SD_raster)
I_coords <- data.frame(ID = 1:nrow(SD_coords), tstep = 1, x_coord = SD_coords[, 1], 
                       y_coord = SD_coords[, 2], progen_ID = NA)

I_coords$secondaries <- rnbinom(nrow(I_coords), mu = 1.2, size = 0.4)
x_min <- min(SD_coords[, 1])
x_max <- max(SD_coords[, 1])
y_max <- max(SD_coords[, 2])
y_min <- min(SD_coords[, 2])
dispersalShape <- 0.3484
dispersalScale <- 41.28/100
t = 1
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
Sys.time()
check <- foreach(i = 1:nrow(I_coords), .errorhandling = 'remove', .combine = rbind
) %dopar% {
  if (I_coords$secondaries[i] > 0) {
    foreach(j = 1:I_coords$secondaries[i], .combine = rbind
            ) %do% { 
      if (j == 1) { # need progenitor coords for 1st movement
        origin_x <- I_coords$x_coord[i]
        origin_y <- I_coords$y_coord[i]
      } else { # progenitors for previous movement
        origin_x <- E_coords$x[nrow(E_coords)]
        origin_y <- E_coords$y[nrow(E_coords)]
      }
      distance <- rweibull(1, shape=dispersalShape, scale=dispersalScale) # Distance to move in km
      angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at
      
      # Convert to a vector and move
      x_new <- (sin(angle) * distance * 1000) + origin_x # convert to m
      y_new <- (cos(angle) * distance * 1000) + origin_y
      
      return(data.frame(tstep = t, x_coord = x_new, y_coord = y_new,  
                                 progen_ID = I_coords$ID[i], path_ID = j))
    }
  }
}
stopCluster(cl)


Sys.time()
for (i in 1:nrow(I_coords)){
  if (I_coords$secondaries[i] > 0) {
    ## foreach look so that we can combine with rbind the E_coords
    for (j in 1:I_coords$secondaries[i]){ 
      if (j == 1) { # need progenitor coords for 1st movement
        origin_x <- I_coords$x_coord[i]
        origin_y <- I_coords$y_coord[i]
      } else {
        origin_x <- E_coords$x[nrow(E_coords)]
        origin_y <- E_coords$y[nrow(E_coords)]
      }
      distance <- rweibull(1, shape=dispersalShape, scale=dispersalScale) # Distance to move in km
      angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at
      
      # Convert to a vector and move
      x_new <- (sin(angle) * distance * 1000) + origin_x # convert to m
      y_new <- (cos(angle) * distance * 1000) + origin_y
      
      E_coords_new <- data.frame(ID = nrow(E_coords) + 1, tstep = t, 
                                 x_coord = x_new, y_coord = y_new,  
                                 progen_ID = I_coords$ID[i], path_ID = j)
      E_coords <- rbind(E_coords, E_coords_new)
    }
  }
}  
Sys.time()
check2 <- rasterize(check[,2:3], SD_raster, fun = function(x,...)length(x))
values(SD_raster) <- runif(ncell(SD_raster), 0, 1)
sus <- SD_raster
values(sus) <- rbinom(length(values(check2)), size = values(check2), prob = values(sus))
## this gives you a raster!

## We need a value from the raster for each point!
check3 <- extract(SD_raster, check[,2:3])
