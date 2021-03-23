# Connectivity metrics ----
conn_stats <- function(names = c("I_dt", "extra_pars", 
                                 "start_date", "days_in_step", 
                                 "V_mat", "N_mat", "tmax",
                                 "cell_ids", "loc_ids",
                                 "prop_start_pop", 
                                 "break_threshold")) {
  
  # Get the objects you need from the environment above this one
  list2env(use_mget(names, envir_num = 2), envir = environment())
  
  I_dt <- I_dt[infected == TRUE]
  
  # first filter all to the minimum tstep (should be passed to extra pars)
  tmin <- extra_pars$tmin
  nlocs <- ncol(V_mat)
  V_mat <- V_mat[, tmin:nlocs]
  N_mat <- N_mat[, tmin:nlocs]
  
  # tabulate I_dt in period of interest
  I_ts <- tabulate(I_dt$t_infectious, nbins = tmax)[tmin:tmax]
  I_ts_local <- tabulate(I_dt[progen_id != -1]$t_infectious, 
                         nbins = tmax)[tmin:tmax]
  I_incs <- I_ts - I_ts_local
  
  # Get consecutive weeks with > 2 * SD above the incursion rate
  nweeks_over_incs <- I_ts > mean(I_incs) + 2 * sd(I_incs)
  if(sum(nweeks_over_incs) > 0) {
    nweeks_over_incs <- max(with(rle(!nweeks_over_incs), lengths[!values]))
  } else {
    nweeks_over_incs <- 0
  }
  
  # Get peak incidence
  peak_inc_total <- max(I_ts)
  peak_inc_local <- max(I_ts_local)
  mean_inc_total <- mean(I_ts)
  mean_inc_local <- mean(I_ts_local)
  
  # Get chain sizes & lengths in time period of interest
  I_gr <- I_dt[progen_id != -1 & t_infectious > tmin][, c("progen_id", "id")]
  
  if(nrow(I_gr) > 0) {
    setnames(I_gr, c("from", "to"))
    I_graph <- graph_from_data_frame(I_gr, directed = TRUE)
    comps <- components(I_graph)
    peak_chain_size <- max(comps$csize)
    mean_chain_size <- mean(comps$csize)
    
    # lengths are a bit trickier
    comp_dt <- data.table(to = as.numeric(names(comps$membership)), 
                          membership = comps$membership)
    I_gr <- comp_dt[I_gr, on = "to"]
    I_gr <- I_gr[, index := paste0(membership, "_", from)]
    chain_lengths <- I_gr[, .(length = uniqueN(index)), by = "membership"]$length
    peak_chain_length <- max(chain_lengths)
    mean_chain_length <- mean(chain_lengths)
    
  } else {
    peak_chain_size <- peak_chain_length <- mean_chain_size <- mean_chain_length <- 0
  }
  
  # village stats
  vill_mat_V <- data.table(loc_ids, V_mat)[, lapply(.SD, sum), by = "loc_ids"]
  vill_mat_N <- data.table(loc_ids, N_mat)[, lapply(.SD, sum), by = "loc_ids"]
  loc_ids <- vill_mat_V$loc_ids
  vill_mat <-   vill_mat_V[, -c("loc_ids")]/vill_mat_N[, -c("loc_ids")]
  conn_met_vill <- apply(vill_mat, 2, conn_metric_vill, 
                         vill_ids = loc_ids, 
                         cov_threshold = extra_pars$cov_threshold * 0.9, 
                         adj_dt = extra_pars$adj_dt)
  max_conn_vill <- max(conn_met_vill)
  mean_conn_vill <- mean(conn_met_vill)
  
  stopped <- prop_start_pop < break_threshold
  
  return(data.table(max_conn_vill, mean_conn_vill, peak_chain_size, 
                    peak_chain_length, mean_chain_size, mean_chain_length,
                    nweeks_over_incs, peak_inc_local, peak_inc_total, 
                    mean_inc_local, mean_inc_total, stopped, 
                    cov_thresh_eff = extra_pars$cov_threshold * 0.9))
  
}

# Function to apply to each column of *vill* level cov matrix ----
conn_metric_vill <- function(x, vill_ids, cov_threshold, adj_dt) {
  
  # filter to cell_ids that have x > threshold
  vills_now <- vill_ids[x >= cov_threshold & x != 0]
  
  # filter adj_dt to these
  adj_now <- adj_dt[from %in% vills_now]
  
  adj_graph <- graph_from_data_frame(adj_now)
  
  npatches <- components(adj_graph)$csize
  conn <- sum(npatches^2)
  return(conn)
  
}

# Function to get neighbor list for sf (only do this once) ----
get_nb_dt <- function(shapefile) {
  
  adj_list <- st_intersects(shapefile)

  rbindlist(lapply(seq_len(length(adj_list)), 
                   function(x) data.table(from = x, to = adj_list[[x]])))
  
}