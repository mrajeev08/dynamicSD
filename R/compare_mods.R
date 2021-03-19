get_reftab_list <- function(dir = "analysis/out/abc_sims") {

  reftabs <- list.files(dir)
  reftabs <- data.frame(fnames = reftabs)
  reftabs %<>%
    tidyr::separate(fnames, into = c("scale", "move", "incs", 
                               "kern", "seed", 
                               NA), "_", 
             remove = FALSE) %>%
    mutate(root_name = paste(scale, move, incs, kern, sep ="_"))
  return(reftabs)
}

read_reftabs <- function(reftab_list, dir = "analysis/out/abc_sims") {
  
  reftl <- mapply(
    function(x, y) {
     out <- fread(paste(dir, x, sep = "/"))
     out$modname <- y
     return(out)
    }, 
    x = reftab_list$fnames, 
    y = reftab_list$root_name,
    SIMPLIFY  = FALSE)

  rbindlist(reftl, fill = TRUE)
  
}

compare_mods <- function(reftable, 
                         par_names = c("R0", "k", "iota"), 
                         exclude = c("stopped", "sim", "break_threshold"),
                         ncores, 
                         paral, 
                         obs_data = NULL, 
                         sampsize = nrow(reftable), 
                         ntree = 500, 
                         predict = TRUE, 
                         return_training = FALSE) {
  
  
  reftab <- reftable[, !c(par_names, exclude), with = FALSE]
  reftab[, modindex := factor(as.numeric(factor(modname)))]
  
  reftab_check <- which(as.vector(reftab[, lapply(.SD, function(x) sum(is.na(x)))]) > 0)
  
  group_list <- as.list(unique(reftab$modindex))
  
  if(length(reftab_check) > 0) {
    reftab <- reftab[, !c(reftab_check), with = FALSE]
  }
  
  if(sampsize < nrow(reftable)) {
    reftab <- reftab[sample(.N, sampsize)]
  } 
  
  modlookup <- data.table(unique.array(reftab[, c("modname", "modindex")]))
  
  reftab <- reftab[, !"modname", with = FALSE]
  
  mod_comp <- abcrf(modindex~., data = reftab, 
                    group = group_list,
                    paral = paral, ncores = ncores, 
                    sampsize = sampsize,
                    ntree = ntree)
  
  if(predict) {
    lda_predicted <- out.abcrf(mod_comp, reftab,
                               predict_ntree = ntree)
  } else {
    lda_predicted <- NULL
  }
  
  if(!is.null(obs_data)) {
    
    mod_predicted <- predict(mod_comp, obs_data, reftab, 
                             paral = paral, ncores = ncores, ntree = ntree)
    
    lda_obs <- data.table(predict(mod_comp$model.lda, obs_data)$x)
    lda_obs[, c("type", 
                "modindex", 
                "post_prob") := .("observed", 
                                  mod_predicted$allocation, 
                                  mod_predicted$post.prob)]
    
    mod_predicted <- data.table(modindex = colnames(mod_predicted$vote), 
                                votes = t(mod_predicted$vote))

  } else {
    mod_predicted <- NULL
  }
  
  err_abcrf <- data.table(err.abcrf(mod_comp, reftab, paral = paral, 
                                    ncores = ncores))
  
  var_imp <- mod_comp$model.rf$variable.importance
  var_imp <- data.table(setNames(stack(var_imp)[2:1], c('summstat','importance')))
  
  if(return_training) {
    training <- reftab
  } else {
    training <- NULL
  }
  return(list(training = training, 
              err = err_abcrf, 
              var_imp = var_imp,
              mod_predicted = mod_predicted, 
              lda_predicted = lda_predicted, 
              lda_obs = lda_obs,
              modlookup = modlookup))
  
}

compare_mod_se <- function(reftable, 
                           par_names = c("R0", "k", "iota"), 
                           exclude = c("stopped", "sim", "break_threshold"),
                           ncores, 
                           paral, 
                           obs_data = NULL, 
                           samp_prop = 0.75, 
                           nsims = 5, 
                           ntree = 500, 
                           predict = TRUE, 
                           return_training = FALSE) {
  
  foreach(i = seq_len(nsims), .combine = comb, 
          .multicombine = TRUE) %do% {
    
    out <- compare_mods(reftable, 
                        par_names, 
                        exclude,
                        ncores, 
                        paral, 
                        obs_data, 
                        sampsize = floor(samp_prop * nrow(reftable)), 
                        ntree, 
                        predict, 
                        return_training)
    invisible(lapply(out, append_col, col_name = "nsim", val = i))
    out[unlist(lapply(out, is.data.table))]
  }
}

