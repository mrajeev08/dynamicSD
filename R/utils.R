# Utility functions for this project --------------------------------------
# MALAVIKA RAJEEV MAY 2016 ------------------------------------------------

# Getting consecutive dates
get.consec <- function (date, format_date, start = "01-01-2002", format_start = "%d-%m-%Y",
                      year1 = 2002, tstep = "bwn", period = FALSE, get.info = FALSE, check = FALSE) {
  
  if (get.info == TRUE) {
    cat ("get.consec
        *For getting consecutive time steps from a set start date* 

        Arguments:
        get.consec (date, format_date, start = '01-01-2002', format_start = '%d-%m-%Y',
                      year1 = 2002, tstep = c('bwn', 'day', 'week', 'month'), period = FALSE, get.info = TRUE, check = FALSE)
        - required: date (your date or vector of dates)
                    format_date (format of your date)
        - defaults: tstep = 'bwn', check = FALSE

        Notes:
        - if tstep = 'bwn', gets date in biweeks, other options are 'day', 'week', 'month'
        - if period = FALSE then uses calendar based functions (i.e. ISOweek, lubridate packages)
          if TRUE, gets consecutive dates by using the interval from the start date (i.e. months = 30.5, biweek = 14, bwn = 7)
        - start =  the start date set to Jan 1st 2002
          format_start = %d-%m-%Y' format of start date input
          year1 = 2002, corresponds to year of start date
        - if check = TRUE, plots distribution of days as a check to make sure function worked

        Dependencies:
        - vars: NA
        - packages: lubridate, ISOweek
        - scripts: none 

        Output:
        - returns: date or vector of dates
        - plots: if check == TRUE, distribution of dates in days as a check to make sure function worked appropriately
         
        ")
  }

  else {

    get.day <- function(date, format_date, start, format_start){
      st <- as.Date (start, format=format_start)
      d <- as.Date(date, format=format_date)
      day.check <- as.numeric(difftime(d, st, units="days"))
    }
    
    day.check <- get.day (date, format_date, start, format_start)
    par(mfrow = c(1, 1))
    par(mar=c(5,5,5,5))
    
    if (check == TRUE) {
      plot (day.check[order(day.check, decreasing=TRUE)], 1:length(day.check), xlab="date check")
    }
      
      if (tstep == "day"){
        return(get.day(date, format_date, start, format_start))
      }
    
      if (tstep == "month"){
      
        if (period == TRUE) {
          get.month.start <- function (date, format_date, start, format_start){
            day <- get.day (date, format_date, start, format_start)
            ceiling (day/30.5)
          }
        return (get.month.start(date, format_date, start, format_start))
      }
    
      else {
        get.month <- function (date, format_date, year1){
          d <- as.Date(date, format=format_date)
          yr <- year(d)-year1
          mon <- month(d) 
          mon + 12*yr
        }
        return (get.month(date, format_date, year1))  
      }
    }

    if (tstep == "bwn"){
      
      if (period == TRUE) {
          get.biweek.start <- function (date, format_date, start, format_start){
            day <- get.day (date, format_date, start, format_start)
            ceiling(day/14)
          }
        return (get.biweek.start(date, format_date, start, format_start))
      }
      
      else {
        get.bwn <- function (date, format_date, year1){
          d <- as.Date(date, format=format_date)
          yr <- year(d)-year1
          wk <- ISOweek (d)
          wk <- (as.numeric(substring (wk, 7, 8)))
          wk[wk==53] <- 52
          wk <- wk + 52*yr
          ceiling(wk/2)
        }
      return(get.bwn(date, format_date, year1))
      }
    }

    if (tstep == "week"){
      if (period == TRUE){
          get.week.start <- function (date, format_date, start, format_start){
            day <- get.day (date, format_date, start, format_start)
            ceiling(day/7)
          }
        return(get.week.start(date, format_date, start, format_start))
      }
      
      else {
        get.week <- function (date, format_date, year1){
          d <- as.Date(date, format=format_date)
          yr <- year(d)-year1
          wk <- ISOweek (d)
          wk <- (as.numeric(substring (wk, 7, 8)))
          wk[wk==53] <- 52
          wk + 52*yr
        }
        return(get.week(date, format_date, year1))
      }
    }
  }
}

# REP VILLNAMES WITH VILLAGE 2002 NAMES -----------------------------------
convert12_02 <- function(data){ # CONVERT VILLAGES FROM 2012 TO 2002! (data 2012)
  data$village <- as.character(data$village) # Convert to character - because otherwise a factor!
  
  # New villages since 2002 census:
  data$village[which(data$village == "Bokore")] <- "Kyambahi" # Bokore from Kyambahi (not Nyichoka)  
  data$village[which(data$village == "Hekwe")] <- "Kenyamonta" # Hekwe from Kenyamonta (not Magatini)
  data$village[which(data$village == "Kenokwe")] <- "Mosongo" # Kenokwe from Mosongo
  data$village[which(data$village == "Kerenero")] <- "Nyamoko" # Kerenero from Nyamoko (not Itununu)
  data$village[which(data$village == "Kitarungu")] <- "Nyansurura" # Kitarungu from Nyansurura
  data$village[which(data$village == "Manyata")] <- "Machochwe" # Manyata from Machochwe
  data$village[which(data$village == "Mbirikiri")] <- "Bonchugu" # Mbirikiri from Bonchugu
  data$village[which(data$village == "Nyamehuru")] <- "Busawe" # Nyamehuru from Busawe
  data$village[which(data$village == "Nyanungu")] <- "Kyambahi" # Nyanungu from Kyambahi (not Nyichoka)
  data$village[which(data$village == "Nyirongo")] <- "Nyamatare" # Nyirongo from Nyamatare
  data$village[which(data$village == "Sogoti")] <- "Kebanchebanche" # Sogoti from Kebanchebanche
  data$village[which(data$village == "Tamkeri")] <- "Mbalibali" # Tamkeri from Mbalibali (not Nyamburi)
  data$village[which(data$village == "Stendi Kuu")] <- "Mugumu" # Stendi Kuu from Mugumu
  # data$village[which(data$village == "Kerukerege")] <- "Maburi" # Kerukerege from Maburi
  # apparenly not, Kerukerege is a kitongoji in Maburi village
  data
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

# Getting probability from rate ---------------------------------------------------------------
get.prob <- function(rate, step) {
  # ## example 1: turn annual waning rate of 0.33 to weekly prob
  # get.prob(rate = 0.33, step = 52)
  # ## example 2: turn annual birth rate of 0.45 to monthly prob
  # get.prob(rate = 0.45, step = 12)
  converted <- (1 + rate)^(1/step) - 1
  return(1 - exp(-converted))
}

