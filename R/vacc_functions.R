##' Simulating vaccination campaign timing at village level
##' ------------------------------------------------------------------------------------------------
#' Simulating vaccination campaign weeks in a given year (1 week/per village)
#' Description
#' Details
#' @param vill_weeks either can generate random week around which each village is vaccinated
#' or assign villages a week before hand passed as vector in same order as vills
#' @return Returned
#' @section Dependencies:
#'     List dependencies here, i.e. packages and other functions
sim.campaigns <- function(data = grid_data, vills = unique(grid_data$villcode),
                          sim_years = 10, burn_in_years = 5, 
                          vill_weeks = sample(1:52, 75, replace = TRUE),...) {
  
  init_week <- sample(1:52, length(vills), replace = TRUE)
  vacc_mat <- matrix(rpois(length(vills)*sim_years, init_week), nrow = length(vills))
  vacc_mat[, 2:ncol(vacc_mat)] <- sapply(2:ncol(vacc_mat), 
                                         function(x) (52*(x - 1) + vacc_mat[,x]))
  
  vacc_weeks <- matrix(1:(sim_years*52), nrow = length(vills), ncol = sim_years*52, byrow = TRUE)
  vaccs <- t(sapply(1:nrow(vacc_weeks), function(x) as.numeric(vacc_weeks[x, ] %in% vacc_mat[x, ])))
  
  ## add in burn in 
  vaccs <- cbind(matrix(0, nrow = length(vills), ncol = burn_in_years*52), vaccs)
  vaccs <- as.data.table(vaccs)
  vaccs$villcode <- vills
  
  vaccs <- as.data.table(list(villcode = data$villcode, 
                                 cell_id = data$cell_id))[vaccs, on = "villcode"]
  vaccs <- vaccs[order(match(vaccs$cell_id, data$cell_id)), ] 
  vaccs <- as.matrix(vaccs[, c("cell_id","villcode"):=NULL]) ## getting rid of id cols
  return(vaccs)
}