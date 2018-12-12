## Profiling
rm(list = ls())

## Libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(rgdal)
library(raster)
library(ISOweek)

## Source in functions
source("functions.R")
source("R/utils.R")

## Data
vacc_data <- read_csv(paste0("data/", list.files("data/", pattern = "Vaccination")))
SD_shape <- readOGR("data/SD_shape/Serengeti_villages UTM_region.shp")
census_data <- read.csv("data/SDcompiled.csv")
pop_data <- read.csv("data/SerengetiPop.csv")

## Get raster data
r <- raster(SD_shape)
res(r) <- 1000
SD_raster <- rasterize(SD_shape, r)

## Data by cellID
grid_data <- get.grid(shapefile = SD_shape, resolution = 1000, pop = pop_data, census = census_data,
                      deaths = 0.44)

## Need for time step functions within data cleaning
y1 <- 2002
ylast <- 2015
steps <- 52
tmax <- (ylast - y1) * steps
start_date <- "01-01-2002"

## Get vacc campaigns at vill level
campaigns <- get.campaigns.WM(vacc = vacc_data, pop = pop_data, shape = SD_shape, threshold = 45)
campaigns$week <- get.consec(campaigns$date_med, format_date = "%Y-%m-%d", start = "01-01-2002", 
                             format_start = "%d-%m-%Y", year1 = 2002, tstep = "week", get.info = FALSE)
ts <- as.data.frame(1:tmax)
names(ts) <- "week"
ts %>% 
  left_join(campaigns) %>%
  dplyr::select(week, total, villcode) %>%
  spread(week, total, fill = 0) %>%
  filter(!is.na(villcode)) -> vacc_mat

# Unhash for testing! -------------------------------------------------------------------------
grid = SD_raster; data = grid_data; vacc = vacc_mat;
I_seed = 2; start_vacc = 0.2;  
R_0 = 1.2; k = 0.2; 
iota = 1; scale_iota = 1;
# weekly number of incursions and scaling factor 
sigma = get.prob(rate = 7/22.3, step = 1); # weekly rate to prob exp -> inf
births = get.prob(rate = grid_data$births, step = 52); # annual birth rate to prob
mu = get.prob(rate = 0.44, step = 52); # annual death rate to prob 
ntimes = tmax; # maximum time step
nu = get.prob(rate = 0.33, step = 52); 
p_revacc = 1; 
# annual waning rate to prob + probability of revaccination
dispersalShape = 0.3484; dispersalScale = 41.28/100;

## No vacc scenario
# vacc_mat[,2:ncol(vacc_mat)] <- ifelse(vacc_mat[,2:ncol(vacc_mat)] > 0, 0, 0)
# start_vacc = 0

# 1. Set-up and tstep 1 -----------------------------------------------------------------------
nlocs <- nrow(data)

## rep incursions (so as to include scaling factor)
inc_rate <- ifelse(length(scale_iota) > 1, iota*scale, rep(iota, length(tmax)))

# state matrices
S <- V <- N <- matrix(NA, nrow = nlocs, ncol = ntimes)
I_dist <- I_all <- E <- matrix(0, nrow = nlocs, ncol = ntimes)

## create coordinates of raster and get min + max values
coords <- coordinates(grid)
x_min <- min(coords[, 1])
x_max <- max(coords[, 1])
y_max <- max(coords[, 2])
y_min <- min(coords[, 2])

## Starting pop + sus 
N[, 1] <- data$start_pop
V[, 1] <- rbinom(n = nlocs, size = data$start_pop, prob = start_vacc)
S[, 1] <- N[, 1] - V[, 1]
sus_rast <- grid


# 2. Simulating infection in first timestep ------------------------------------------------------
# Seeding cases (as if from outside!)
# pick random location to start I seeds
I_dist[ ,1] <- 0
I_locs <- sample((1:nlocs)[which(!is.na(data$start_pop))], I_seed)
I_dist[I_locs, 1] <- 1

## Create infected coordinate data frame to grow
I_coords <- data.frame(ID = 1:I_seed, tstep = rep(1, I_seed), x_coord = coords[I_locs, 1], 
                       y_coord = coords[I_locs, 2], progen_ID = rep(0, I_seed),
                       path_ID = rep(NA, I_seed), sus = NA, trans = NA,
                       infectious = rep(1, I_seed), secondaries = NA)

## Create exposed coordinate data frame to grow
E_coords <- data.frame(ID = NA, tstep = NA, x_coord = NA, y_coord = NA, 
                       progen_ID = NA, path_ID = NA, sus = NA, trans = NA, 
                       infectious = NA, secondaries = NA)
out <- E_coords ## for keeping track of ones that go outside and tstep E_coords

I_coords$secondaries <- rnbinom(nrow(I_coords), mu = R_0, size = k)
counter <- length(I_seed) ## for keeping track of ids

for (i in 1:nrow(I_coords)){
  if (I_coords$secondaries[i] > 0) {
    for (j in 1:I_coords$secondaries[i]){ 
      counter <- counter + 1
      if (j == 1) { # need progenitor coords for 1st movement
        origin_x <- I_coords$x_coord[i]
        origin_y <- I_coords$y_coord[i]
      } else {
        origin_x <- E_coords$x[nrow(E_coords)]
        origin_y <- E_coords$y[nrow(E_coords)]
      }
      distance <- rweibull(1, shape = dispersalShape, 
                           scale = dispersalScale) # Distance to move in km
      angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at
      
      # Convert to a vector and move
      x_new <- (sin(angle) * distance * 1000) + origin_x # convert to m
      y_new <- (cos(angle) * distance * 1000) + origin_y
      
      E_coords_new <- data.frame(ID = counter, tstep = 1, 
                                 x_coord = x_new, y_coord = y_new,  
                                 progen_ID = I_coords$ID[i], path_ID = j,
                                 sus = NA, trans = NA, infectious = NA, secondaries = NA)
      E_coords <- rbind(E_coords, E_coords_new)
    } 
      
    ## Potentially don't need this as we wont have sus for the points outside the district
    out <- rbind(out, filter(E_coords, x_coord <= x_min | x_coord >= x_max |
                               y_coord <= y_min | y_coord >= y_max))
    E_coords_now <- filter(E_coords, tstep == 1, x_coord >= x_min, x_coord <= x_max,
                           y_coord >= y_min, y_coord <= y_max) ## filter out ones outside
    values(sus_rast) <- S[, 1]/N[, 1] ## probability that contact will be with a susceptible
    
    E_coords_now$sus <- extract(sus_rast, cbind(E_coords_now$x_coord, E_coords_now$y_coord))
    E_coords_now$trans <- rbinom(nrow(E_coords_now), size = 1, prob = E_coords_now$sus)
    E_coords_now %>%
      drop_na(trans) %>%
      filter(trans == 1) -> E_coords
  } else
    E_coords_now <- filter(E_coords, tstep == -1)
}

if(nrow(E_coords) > 1) {
  E[, 1] <- values(rasterize(cbind(E_coords$x_coord, E_coords$y_coord), grid, 
                             fun= function(x,...)length(x)))
  E[is.na(E[,1]), 1] <- 0
} else {
  E[, 1] <- 0
}

## Finish subtracting everyone out
S[, 1] <- S[, 1] - E[ ,1] ## subtract out newly exposed guys from S

# 3. Simulate for rest of time steps -------------------------------------------------------------
for (t in 2:ntimes) {
  ## t = 2
  print(paste(t, "/", tmax, "weeks"))
  # 3a. Simulate vaccination ------------------------------------------------------------------------
  ## TO DO: ## speed up not using filter
  ## TO DO: ## make sure leftover vaccinated goes to adjacent vills
  data %>%
    mutate(S = S[, t-1], V = V[, t-1]) %>%
    dplyr::select(villcode, S, V) %>%
    group_by(villcode) %>%
    summarize_all(sum, na.rm = TRUE) -> now
  vacc %>%
    dplyr::select(villcode, vacc = t + 1) %>%
    left_join(now) %>%
    mutate(p_vacc = ifelse(vacc - p_revacc*V <= 0, vacc/S,
                           (vacc - p_revacc*V)/S)) %>%
    mutate(p_vacc = ifelse(p_vacc > 1, 1, p_vacc)) %>% 
    right_join(data) %>%
    mutate(n_vacc = rbinom(nlocs, size = S, prob = p_vacc*prop_pop)) %>%
    replace_na(list(n_vacc = 0)) -> vacc_now
  
  # 3b. Sequential competing transitions for S + V -----------------------------------------------
  ## S transitions, sequenctially for competing probs (-vacc, - deaths)
  S[, t] <- S[, t-1] - vacc_now$n_vacc ## - vaccinated
  S[, t] <- S[, t] - rbinom(nlocs, size = S[, t], prob = mu)
  
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
  
  E_coords %>% 
    drop_na(infectious) %>%
    filter(infectious == 1) -> I_coords_in ## add any that became infectious to I_coords_now
  
  ## Only count exposed -> infectious from within to subtract out
  if(nrow(I_coords_in) > 0) {
    I_dist[, t] <- values(rasterize(cbind(I_coords_in$x_coord, I_coords_in$y_coord), grid, 
                             fun= function(x,...)length(x)))
    I_dist[is.na(I_dist[, t]), t] <- 0 
  } else {
    I_dist[, t] <- 0
  }
  
  E_coords %>% 
    drop_na(infectious) %>%
    filter(infectious == 0) -> E_coords ## remove any that became infectious from E_coords
  
  # 3d. Add in incursions --------------------------------------------------------------------
  ## incursions (To do? : could add in option to weight by distance from edge)
  incursions <- rpois(1, inc_rate) # number of incursions in week
  if (incursions > 0){
    start_counter <- counter + 1 
    counter <- counter + incursions
    I_all[, t] <- 0
    
    ## only happen in populated places
    I_locs <- sample((1:nlocs)[which(!is.na(data$start_pop))], incursions) 
    I_all[I_locs, t] <- 1 
    
    I_coords_out <- data.frame(ID = start_counter:counter, tstep = t, x_coord = coords[I_locs, 1],
                                y_coord = coords[I_locs, 2], progen_ID = 0, path_ID = NA, sus = NA, 
                                trans = NA, infectious = 1, 
                                secondaries = NA)
    I_coords_now <- rbind(I_coords_in, I_coords_out)
  }
  
  I_all[, t] <-I_all[, t] + I_dist[, t]
  
  # 3e. Simulate individual based contacts + movement --------------------------------------------
  E_coords_now <- filter(E_coords, tstep == t) ## to get empty df to bind to!
  
  if(nrow(I_coords_now) > 0){
    
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
          check <- NA
          while (is.na(check)) {
            distance <- rweibull(1, shape = dispersalShape, 
                               scale = dispersalScale) # Distance to move in km
            angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at
          
            # Convert to a vector and move
            x_new <- (sin(angle) * distance * 1000) + origin_x # convert to m
            y_new <- (cos(angle) * distance * 1000) + origin_y
            
            ## This means will only go somewhere within a vill
            check <- extract(grid, cbind(x_new, y_new))
            ## TO DO ##: Need some code in here so that it goes to populated place!
          }
          
          E_coords_new <- data.frame(ID = counter, tstep = t, 
                                     x_coord = x_new, y_coord = y_new,  
                                     progen_ID = I_coords_now$ID[i], path_ID = j, 
                                     sus = NA, trans = NA, 
                                     infectious = NA, secondaries = NA)
          last_coords <- c(x_new, y_new)
          E_coords_now <- rbind(E_coords_now, E_coords_new)
        }
      }
    }
  } else {
    E_coords_now <- I_coords_now ## gets you empty df
  }
  
  # out <- rbind(out, filter(E_coords_now, x_coord <= x_min | x_coord >= x_max |
  #                            y_coord <= y_min | y_coord >= y_max)) # store ones that went outside
  # E_coords_now <- filter(E_coords_now, x_coord >= x_min, x_coord <= x_max,
  #                   y_coord >= y_min, y_coord <= y_max) ## filter out ones where S[,t] is NA
  
  # 3f. Aggregate probability that contact is with susceptible -----------------------------------
  
  if(nrow(E_coords_now) > 0) {
    values(sus_rast) <- S[, t]/(S[, t] + E[, t-1] + V[, t])   
    ## probability that contact will be with a susceptible, E[,t-1] = prob with I_dist[,t] too
    E_coords_now$sus <- extract(sus_rast, cbind(E_coords_now$x_coord, E_coords_now$y_coord))
    E_coords_now$trans <- rbinom(nrow(E_coords_now), size = 1, prob = E_coords_now$sus)
    E_coords_now %>%
      drop_na(trans) %>%
      filter(trans == 1) -> E_coords_now
  } else {
    E_coords_now <- E_coords_now # to get empty df to bind!
  }
  
  if (nrow(E_coords_now) > 0) {
    E_new <- values(rasterize(cbind(E_coords_now$x_coord, E_coords_now$y_coord), 
                              grid, fun= function(x,...)length(x)))
    E_new[is.na(E_new)] <- 0
  } else {
    E_new <- rep(0, nlocs)
  }
  
  ## Finish subtracting everyone out
  E[, t] <- E[, t-1] + E_new - I_dist[, t] ## - infectious from within district + newly exposed
  S[, t] <- S[, t] - E_new ## subtract out newly exposed guys in the week
  N[, t] <- S[, t] + V[, t] + E[, t] + I_dist[, t]
  
  E_coords <- rbind(E_coords, E_coords_now)
  
  I_coords <- rbind(I_coords, I_coords_now) ## to build trees!
}

sum_times <- function(vector, steps, na.rm=TRUE) {    # 'matrix'
  nv <- length(vector)
  if (nv %% steps)
    vector[ceiling(nv / steps) * steps] <- NA
  colSums(matrix(vector, steps), na.rm = na.rm)
}

mIobs <- apply(I_all[,1:(ncol(I_all)-3)], 1 , sum_times, steps = 4)
plot(rowSums(mIobs[ ,!is.na(grid_data$village_ID)], na.rm = TRUE), type = "l")
plot(colSums(I_all[!is.na(grid_data$village_ID), ], na.rm = TRUE), type = "l")

