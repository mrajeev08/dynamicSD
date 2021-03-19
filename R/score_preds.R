# scoring rules & summary functions ----
get_tempstats <- function(out, obs_data, quants = c(0.5, 0.9), 
                          nsamp = 100, ncurves = 100, nbest = 5) {
  
  out_sims <- copy(out)
  
  preds <- as.matrix(dcast(out_sims, cal_month ~ sim, value.var = "I_ts")[, -1])
  if(nbest > max(out_sims$sim)) nbest <- max(out_sims$sim)
     
  stats <- as.data.table(lapply(list(min, 
                                     max), 
                                function(x) apply(preds, 1, x)))
  setnames(stats, c("min", "max"))
  
  curve_scores <- get_centrality(preds_dt = out_sims, 
                                 nsamp = nsamp, 
                                 ncurves = ncurves)
  
  best_sims <- mapply(get_repsims, 
                      preds_dt = list(out_sims),
                      curve_scores = list(curve_scores),
                      type = c("all", "times"),
                      n = nbest, 
                      SIMPLIFY = FALSE)

  envelopes <- mapply(get_envelope, 
                      type = c("all", "times"),
                      quantile = quants,
                      preds_dt = list(out_sims),
                      curve_scores = list(curve_scores), 
                      SIMPLIFY = FALSE)
  
  out_dts <- Reduce(function(...) merge(..., all = TRUE), 
                    c(envelopes, best_sims))
  if(!is.null(obs_data)) {
    
    scores <- as.data.table(lapply(list(crps_sample, 
                                        logs_sample, 
                                        dss_sample), 
                                   function(x) x(obs_data$cases_by_month, 
                                                 dat = preds)))
    setnames(scores, c("crps_score", "logs_score", "dss_score"))
    
    # also get sim closest to data
    scores[, paste0("best_dat_", 1:nbest) := rank_sims(obs_data, 
                                                        preds, 
                                                        n = nbest)]
    
  } else {
    scores <- NULL
  }
  
  cbind(out_dts, stats, scores)

}

get_repsims <- function(preds_dt, curve_scores, 
                         n = 5,
                         type = c("all", "times")) {
  type <- match.arg(type)
  curve_scores[, score := get(paste("score", type, sep = "_"))]
  setorder(curve_scores, score)
  sims_in <- curve_scores$sim[1:n]
  sims_out <- dcast(preds_dt[sim %in% sims_in], 
                        cal_month ~ sim, value.var = "I_ts")
  setnames(sims_out, 2:(n + 1), paste("best", type, 1:n, sep = "_"))
  
}

get_envelope <- function(preds_dt, curve_scores, quantile= 0.75, 
                         type = c("all", "times")) {
  
  # sort sims by scores
  type <- match.arg(type)
  curve_scores[, score := get(paste("score", type, sep = "_"))]
  setorder(curve_scores, score)
  thresh_ind <- round(max(curve_scores$sim)*quantile)
  sims_in <- curve_scores$sim[1:thresh_ind]
  env_all <- preds_dt[sim %in% sims_in][, .(min = min(I_ts),
                                             max = max(I_ts)), 
                                        by = cal_month]
  setnames(env_all, c("min", "max"), 
           paste(c("min", "max"), type, quantile, sep = "_"))

  
  return(env_all)
  
}

get_centrality <- function(preds_dt, nsamp = 100, ncurves = 100) {
  
  out <- 
    foreach(i = seq_len(nsamp), .combine = 'rbind') %dopar% {
    
    smp_inds <- sample(max(preds_dt$sim), ncurves, 
                       replace = ifelse(max(preds_dt$sim) < ncurves, 
                                        TRUE, 
                                        FALSE))
    
    preds_now <- preds_dt[sim %in% smp_inds][, .(min = min(I_ts),
                                                 max = max(I_ts)), 
                                             by = cal_month]
    pred_scores <- preds_dt[preds_now, on = "cal_month"]
    
    pred_scores[, .(within = sum(I_ts <= max & I_ts >= min), 
                    all = all(I_ts <= max & I_ts >= min)), by = sim]
  }
  
  out[, .(score_all = sum(within), score_times = sum(all)), by = sim]

  
}

rank_sims <- function(obs_data, preds, n) {
  
  rank_rmse <- colSums(sqrt(apply(preds, 2, 
                                  function(x) x - obs_data$cases_by_month)^2))/nrow(preds)
  
  data.table(preds[, order(rank_rmse)[1:n]])

}
