
# Getting probability from rate ---------------------------------------------------------------
get.prob <- function(rate, step) {
  # ## example 1: turn annual waning rate of 0.33 to weekly prob
  # get.prob(rate = 0.33, step = 52)
  # ## example 2: turn annual birth rate of 0.45 to monthly prob
  # get.prob(rate = 0.45, step = 12)
  converted <- (1 + rate)^(1/step) - 1
  return(1 - exp(-converted))
}

# Aggregating village vaccination to campaigns -----------------------------------------------
## To do: need to add in part to check ones with no Vill 2002 names!
## and shout warning of how many doggies we're losing to this
get.campaigns.WM <- function(vacc = vacc_data, pop = pop_data, threshold, shape) {
  ## vacc = vacc_data; pop = pop_data; shape = SD_shape; threshold = 45;
  vacc <- filter(vacc, District == "Serengeti")
  vacc$villcode <- as.factor(pop$VILLCODES[match(vacc$`Village 2002`, pop$Village_2002)])
  
  cat("Step 1: fixing by point locations:\n")
  print(vacc$Village[is.na(vacc$villcode)])
  loc_match <- filter(vacc, is.na(villcode))
  loc_match$villcode <- extract(shape, cbind(loc_match$`UTM Easting`, loc_match$`UTM Northing`))$VILLCODE
  vacc$villcode[is.na(vacc$villcode)] <- loc_match$villcode[match(vacc$Reference[is.na(vacc$villcode)], 
                                                                  loc_match$Reference)]
  
  cat("Step 2: fixing by fuzzy matching vill 2012 to vill 2002 names:\n")
  print(vacc$Village[is.na(vacc$villcode)])
  
  check <- unique(vacc$Village[is.na(vacc$villcode)])
  for (i in 1:length(check)){
    vacc$villcode[vacc$Village == check] <- pop$VILLCODES[agrep(check[i], pop$Village_2012)[1]]
  }
  
  cat("Warning: if greater than 0, then still unmatched vaccinations:\n")
  print(vacc$Village[is.na(vacc$villcode)])
  
  vacc %>%
    mutate(cats = replace_na(`Cats Vaccinated`, 0), # merge up doses + dogs + cats for total 
           dogs = replace_na(`Dogs Vaccinated`, 0), 
           total = coalesce(Doses, (dogs + cats)), 
           date = dmy(`Actvitity Date`)) %>%  # format date 
    group_by(villcode, date) %>% # ordered by date
    summarize(total = sum(total)) %>%  # get sum if same date in same village
    spread(villcode, total) %>%
    mutate(diff = c(diff(date), diff(date)[length(date)-1]),
           group = if_else(diff > threshold, 1, 0),
           group = rep(1:length(rle(group)$values), rle(group)$lengths)) %>%
    dplyr::select(-diff) %>% ## need to look at where these are and assign villages!
    gather(villcode, doses, -group, -date) %>%
    group_by(group, villcode) %>%
    filter(!is.na(doses)) %>%
    summarize(total = sum(doses), 
              date_med = min(date) + (max(date) - min(date))/2) -> campaigns 
  return(campaigns)
}

# Getting gridded data -----------------------------------------------------
## Gets the village ID, starting pop, the proportion of the pop by village, 
## HDR, and the birth rate for each cell in grid
get.grid <- function(shapefile, resolution = 1000, pop = pop_data, census = census_data, 
                     deaths = 0.44) {
  # pop = pop_data; census = census_data; deaths = 0.44; shapefile = SD_shape;
  
  ## 1. Get raster
  r <- raster(shapefile)
  res(r) <- 1000
  SD_raster <- rasterize(shapefile, r)
  rast_attributes <- SD_raster@data@attributes[[1]]
  shapefile$ID_match <-rast_attributes$ID[match(shapefile$VILLCODE, rast_attributes$VILLCODE)]
  
  ## 2. Get proportion of dogs in each cell
  census_dogs <- rasterize(cbind(census$X, census$Y), r,
                           field = census$dogs + census$pups, fun = sum)
  summarize_dogs <- as.data.frame(cbind(census_dogs@data@values, SD_raster@data@values))
  summarize_dogs %>%
    group_by(village_ID = V2) %>%
    summarize(dogs_total_census = sum(V1, na.rm = TRUE)) %>%
    right_join(shapefile@data, by = c("village_ID" = "ID_match")) -> shapefile@data
  raster_pop <- rasterize(shapefile, r, field = shapefile$dogs_total_census)
  prop <- census_dogs/raster_pop
  
  # 3. Get starting pop and births from village censuses 
  pop$VILLCODES <- as.factor(pop$VILLCODES)
  shapefile@data %>%
    left_join(pop, by = c("VILLCODE" = "VILLCODES")) %>%
    group_by(Village_2002) %>%
    summarize(pop2002 = sum(Population.2002, na.rm = TRUE), 
              pop2012 = sum(Population.2012, na.rm = TRUE), 
              villcode = VILLCODE[1]) %>%
    mutate(growth = exp(log(pop2012/pop2002)/10)) %>%
    right_join(shapefile@data, by = c("villcode" = "VILLCODE")) -> shapefile@data
  shapefile$births <- deaths + (shapefile$growth - 1)

  ## 4. Get HDR and starting dog population from starting vill population
  census_humans <- rasterize(cbind(census$X, census$Y), r,
                             field = census$Adults + census$Children, fun = sum)
  HDR <- census_humans/census_dogs
  births <- rasterize(shapefile, r, field = shapefile$births)
  start_pop <- rasterize(shapefile, r, field = shapefile$pop2002)
  values(start_pop) <- rbinom(ncell(start_pop), size = values(start_pop), 
                              prob = 1/values(HDR)*values(prop))
  
  ## 5. Return stacked raster file
  dat <- as_tibble(cbind(SD_raster@data@values, HDR@data@values, births@data@values, 
                             prop@data@values, start_pop@data@values))
  colnames(dat) <- c("village_ID", "HDR", "births", "prop_pop", "start_pop")
  dat$villcode <- rast_attributes$VILLCODE[match(dat$village_ID, rast_attributes$ID)]
  dat$cell_id <- 1:nrow(dat) ## in order to use identity of cell!
  return(dat)
}


# Simulating with multinomial probability (not summing to 1) ----------------------------------
## Needs to be a matrix!
sim.multinom <- function(sizes = c(10, 15, 20, 30), 
                         probs = cbind(0.05, 0.01)) {
  trans <- matrix(NA, nrow = length(sizes), ncol = ncol(probs) + 1)
  
  if (nrow(probs) == 1){
    prob_vec <- c(probs[1, ], 1 - sum(probs[1, ]))
    for (i in 1:length(sizes)){
      if(!is.na(sizes[i])) {
        trans[i, ] <- rmultinom (1, size = sizes[i], prob_vec)
      }
    }
  }
  
  if(nrow(probs) > 1) {
    for (i in 1:length(sizes)){
      if(!is.na(sizes[i])) {
        prob_vec <- c(probs[i, ], 1 - sum(probs[i, ]))
        trans[i, ] <- rmultinom (1, size = sizes[i], prob_vec)
      }
    }
  }
  return(trans) 
}

# Growing data frames -----------------------------------------------------------------------
## here elems is a named list!
## only at super high # of rows does this become more efficient
create.df <- function(elems){
  dt <- as.data.table(elems)
  setattr(dt, 'rowcount', max(lengths(elems)))
}

append.df <- function(dt, elems) {
  n <- attr(dt, 'rowcount')
  
  if(is.null(n)) {
    n <- nrow(dt)
  }
  
  max_row <- n + max(lengths(elems))
  
  if(max_row > nrow(dt)) {
    tmp <- elems[(nrow(dt)+1):(nrow(dt)*2)]
    dt <- rbindlist(list(dt, tmp), fill = TRUE, use.names = TRUE)
  }
  
  dt[(n+1):max_row, ] <- elems
  
  setattr(dt, 'rowcount', max_row)

  return(dt)
}

access.df <- function(dt) {
  n <- attr(dt, 'rowcount')
  return(as.data.table(dt[1:n, ]))
}

