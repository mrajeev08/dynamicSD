estimate_pars <- function(reftable, 
                          par_names = c("R0", "k", "iota"), 
                          exclude = c("stopped", "sim", "break_threshold"),
                          ncores, 
                          paral, 
                          obs_data = NULL, 
                          sampsize = nrow(reftable), 
                          ntree = 500, 
                          predict = TRUE, 
                          predict_nsimul = 0,
                          return_training = FALSE) {
  
    
    reftab <- reftable[, !c(exclude), with = FALSE]
    
    modname <- unique(reftab$modname)
    
    if(length(modname) > 1) stop("Attempting to estimate params with more than one model!")

    reftab <- reftab[, !"modname", with = FALSE]
    
    if(sampsize < nrow(reftable)) {
      reftab <- reftab[sample(.N, sampsize), ]
    } 
      
    out <- 
      foreach(i = seq_len(length(par_names)), .combine = comb, 
            .multicombine = TRUE) %do% {
        
        reftab_i <- reftab[, !par_names[-i], with = FALSE]
        setnames(reftab_i, par_names[i], "par")
        
        if(predict_nsimul > 0) {
          inds <- 1:predict_nsimul
          simuls <- reftab_i[inds, ]
          reftab_i <- reftab_i[-inds, ]
        }
        
        par_obj <- regAbcrf(par~., data = reftab_i, ntree = ntree, 
                            paral = paral, ncores = ncores)
        
        if(predict & !is.null(obs_data)) {

          out_preds <- out.regAbcrf(par_obj, obs = obs_data, 
                                    training = reftab_i, 
                                    paral = paral, 
                                    ncores = ncores)
          out_preds[, c("sim_id", "param") := .(1:.N, par_names[i])]
          
          out_stats <- predict(par_obj, obs = obs_data, 
                               training = reftab_i, 
                               paral = paral, ncores = ncores,
                               ntree = ntree)
          
          out_stats <- as.data.table(out_stats)
          setnames(out_stats, c("expectation", "median", "variance_postmse", 
                                "variance_cdf", "quantile_0.025", 
                                "quantile_0.975", "post_nmae_mean"))
         
        } else {
          out_preds <- out_stats <- NULL
        }
        
        if(predict_nsimul > 0) {
          preds_simul <- predict(par_obj, obs = simuls, 
                               training = reftab_i, 
                               paral = paral, ncores = ncores,
                               ntree = ntree)
          preds_simul <- as.data.table(preds_simul)
          setnames(preds_simul, c("expectation", "median", "variance", 
                                "variance_cdf", "quantile_0.025", 
                                "quantile_0.975", "post_nmae"))
          preds_simul$param <- par_names[i]
          preds_simul$true_val <- simuls$par
        } else {
          preds_simul <- NULL
        }
        err_abcrf <- data.table(err.regAbcrf(par_obj, reftab_i, paral = paral, 
                                             ncores = ncores))
        
        var_imp <- par_obj$model.rf$variable.importance
        var_imp <- data.table(setNames(stack(var_imp)[2:1], 
                                       c('summstat','importance')))
        out <- list(preds = out_preds, err = err_abcrf, var_imp = var_imp, 
                    stats = out_stats, preds_test = preds_simul)
        
        out <- lapply(out, append_col, val = par_names[i], col_name = "param")
        out
      }
    
    out <- lapply(out, append_col, val = modname, col_name = "modname")
    
    if(return_training) {
      return(list(out = out, training = training))
    } else {
      return(out)
    }
    
}

estimate_par_se <- function(reftable, 
                           par_names = c("R0", "k", "iota"), 
                           exclude = c("stopped", "sim", "break_threshold"),
                           ncores, 
                           paral, 
                           obs_data = NULL, 
                           ntree = 500, 
                           predict = TRUE, 
                           samp_prop = 0.75, 
                           nsims = 5,
                           predict_nsimul = 0) {
  
  foreach(i = seq_len(nsims), .combine = comb, 
          .multicombine = TRUE) %do% {
            
            out <- estimate_pars(reftable, 
                                  par_names,
                                  exclude,
                                  ncores, 
                                  paral, 
                                  obs_data, 
                                  sampsize = floor(samp_prop * nrow(reftable)), 
                                  ntree, 
                                  predict, 
                                  predict_nsimul,
                                  return_training = FALSE)
            
            out <- lapply(out, append_col, col_name = "nsim", val = i)
            out
  } -> out
  
  return(out)
  
}

