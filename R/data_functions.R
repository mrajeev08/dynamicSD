# Aggregating village vaccination to campaigns -----------------------------------------------
## To do: need to add in part to check ones with no Vill 2002 names!
## and shout warning of how many doggies we're losing to this
## and add adjacency correction for villages!
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

# Getting demographic data -----------------------------------------------------
## Gets the village ID, starting pop, the proportion of the pop by village, 
## HDR, and the growth rate for each cell in grid
## Gets the village growth rate and starting pop for each village

get.demdata <- function(shapefile, res_meters = 1000, pop = pop_data, census = census_data) {
  ## To test code
  pop = pop_data; census = census_data; res_meters = 1000; deaths = 0.44; shapefile = SD_shape;
  
  # 1. Get starting pop and births from village censuses 
  pop$VILLCODES <- as.factor(pop$VILLCODES)
  shapefile@data %>%
    left_join(pop, by = c("VILLCODE" = "VILLCODES")) %>%
    group_by(Village_2002) %>%
    summarize(pop2002 = sum(Population.2002, na.rm = TRUE), 
              pop2012 = sum(Population.2012, na.rm = TRUE), 
              villcode = VILLCODE[1]) %>%
    mutate(growth = exp(log(pop2012/pop2002)/10)) %>%
    right_join(shapefile@data, by = c("villcode" = "VILLCODE")) -> shapefile@data
  
  if (is.numeric(res_meters)) {
    ## 2. Get raster
    r <- raster(shapefile)
    res(r) <- res_meters
    SD_raster <- rasterize(shapefile, r)
    rast_attributes <- SD_raster@data@attributes[[1]]
    shapefile$ID_match <-rast_attributes$ID[match(shapefile$villcode, rast_attributes$villcode)]
    
    ## 3a. Get proportion of dogs in each cell (from census)
    census_dogs <- rasterize(cbind(census$X, census$Y), r,
                             field = census$dogs + census$pups, fun = sum)
    summarize_dogs <- as.data.frame(cbind(census_dogs@data@values, SD_raster@data@values))
    summarize_dogs %>%
      group_by(village_ID = V2) %>%
      summarize(dogs_total = sum(V1, na.rm = TRUE)) %>%
      right_join(shapefile@data, by = c("village_ID" = "ID_match")) -> shapefile@data
    raster_pop <- rasterize(shapefile, r, field = shapefile$dogs_total)
    prop <- census_dogs/raster_pop
    
    ## 3b.And also dogs vaccinated
    census_vacc <- rasterize(cbind(census$X, census$Y), r,
                             field = census$vacc_dogs + census$vacc_pups, fun = sum)
    cov <- census_vacc/census_dogs
    ##' add small prob so that allocates equally if vaccinations happen in place where
    ##' no coverage reported in census and can be allocated equally when 'leftovers'
    cov[cov > 1] <- 1 ## max of 1 
    cov[is.na(cov) | cov == Inf] <- 0
    cov <- cov + 1e-5 
    
    ## 4. Get HDR and starting dog population from starting vill population
    census_humans <- rasterize(cbind(census$X, census$Y), r,
                               field = census$Adults + census$Children, fun = sum)
    summarize_humans <- as.data.frame(cbind(census_humans@data@values, SD_raster@data@values))
    summarize_humans %>%
      group_by(village_ID = V2) %>%
      summarize(humans_total = sum(V1, na.rm = TRUE)) %>%
      right_join(shapefile@data, by = c("village_ID" = "village_ID")) -> shapefile@data
    shapefile$HDR <- shapefile$humans_total/shapefile$dogs_total
    shapefile$start_pop <- shapefile$pop2002/shapefile$HDR
    HDR <- rasterize(shapefile, r, field = shapefile$HDR)
    growth <- rasterize(shapefile, r, field = shapefile$growth) 
    start_pop <- rasterize(shapefile, r, field = shapefile$start_pop) ## just the starting vill pops
    
    ## 5. Return stacked raster file
    dat <- as.data.table(list(village_ID = SD_raster@data@values, 
                              HDR = HDR@data@values, growth = growth@data@values, 
                              prop_pop = prop@data@values, start_pop = start_pop@data@values, 
                              cov = cov@data@values))
    dat$start_pop <- rbinom(nrow(dat), size = round(dat$start_pop), prob = dat$prop_pop)
    dat$villcode <- rast_attributes$villcode[match(dat$village_ID, rast_attributes$ID)]
    dat$cell_id <- 1:nrow(dat) ## in order to use identity of cell!
    return(dat)
  } else {
    ## For village data get HDR + start_pop
    census %>%
      group_by(Village.2002) %>%
      summarize(dogs_total = sum(dogs + pups, na.rm = TRUE), 
                humans_total = sum(Adults + Children, na.rm = TRUE), 
                HDR = humans_total/dogs_total) %>%
      right_join(shapefile@data, by = c("Village.2002" = "Village_2002")) %>% # keep all shapefile data
      mutate(start_pop = pop2002/HDR) -> shapefile@data
    return(shapefile@data)
  }
}

# Cleaning rabid  -----------------------------------------------------------------------------


