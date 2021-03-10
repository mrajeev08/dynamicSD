sample_posterior <- function(N = 1000, 
                             object, obs, training, paral = FALSE, 
                             ncores= if(paral) max(detectCores()-1, 1) else 1) {
 
      # Checking args
      if (!inherits(object, "regAbcrf")) 
        stop("object not of class regAbcrf")
      
      if (!inherits(training, "data.frame"))
        stop("training needs to be a data.frame object")
      
      if (!inherits(obs, "data.frame")) 
        stop("obs needs to be a data.frame object")
  
      if (nrow(obs) == 0L || is.null(nrow(obs)))
        stop("no data in obs")
  
      if (nrow(training) == 0L || is.null(nrow(training)))
        stop("no simulation in the training reference table (response, sumstat)")
  
      if ( (!is.logical(paral)) || (length(paral) != 1L) )
        stop("paral should be TRUE or FALSE")
  
      if(is.na(ncores)){
        warning("Unable to automatically detect the number of CPU cores, 
                \n1 CPU core will be used or please specify ncores.")
        ncores <- 1
      }
  
      x <- obs
      
      if(!is.null(x)){
        if(is.vector(x)){
          x <- matrix(x, ncol = 1)
        }
        if (nrow(x) == 0) 
          stop("obs has 0 rows")
        if (any(is.na(x))) 
          stop("missing values in obs")
      }
      
      # resp and sumsta recover
      mf <- match.call(expand.dots=FALSE)
      mf <- mf[1]
      mf$formula <- object$formula
      
      mf$data <- training
      
      training <- mf$data
      
      mf[[1L]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame() )
      mt <- attr(mf, "terms")
      resp <- model.response(mf)
      
      obj <- object$model.rf
      inbag <- matrix(unlist(obj$inbag.counts, use.names=FALSE), ncol=obj$num.trees, byrow=FALSE)
      
      obj[["origNodes"]] <- predict(obj, training, predict.all=TRUE, num.threads=ncores)$predictions
      obj[["origObs"]] <- model.response(mf)
      
      #####################
      
      origObs <- obj$origObs
      origNodes <- obj$origNodes
      
      nodes <- predict(obj, x, predict.all=TRUE, num.threads=ncores )$predictions
      if(is.null(dim(nodes))) nodes <- matrix(nodes, nrow=1)
      ntree <- obj$num.trees
      nobs <- object$model.rf$num.samples
      nnew <- nrow(x)
      weights <- abcrf:::findweights(origNodes, nodes, inbag, 
                                     as.integer(nobs), 
                                     as.integer(nnew), 
                                     as.integer(ntree)) # cpp function call
      
      weights.std <- weights/ntree
      
      post_sample <- matrix(NA, ncol = nnew, nrow = N)
      
      for(i in seq_len(nnew)) {
        post_sample[, i] <- sample(resp, N, prob = weights.std[, i], 
                                   replace = TRUE)
      }
      
      return(post_sample)
      
  }

posterior_density <- function(object, obs, training, paral = FALSE, 
                             ncores= if(paral) max(detectCores()-1, 1) else 1) {
  
  # Checking args
  if (!inherits(object, "regAbcrf")) 
    stop("object not of class regAbcrf")
  
  if (!inherits(training, "data.frame"))
    stop("training needs to be a data.frame object")
  
  if (!inherits(obs, "data.frame")) 
    stop("obs needs to be a data.frame object")
  
  if (nrow(obs) == 0L || is.null(nrow(obs)))
    stop("no data in obs")
  
  if (nrow(training) == 0L || is.null(nrow(training)))
    stop("no simulation in the training reference table (response, sumstat)")
  
  if ( (!is.logical(paral)) || (length(paral) != 1L) )
    stop("paral should be TRUE or FALSE")
  
  if(is.na(ncores)){
    warning("Unable to automatically detect the number of CPU cores, 
                \n1 CPU core will be used or please specify ncores.")
    ncores <- 1
  }
  
  x <- obs
  
  if(!is.null(x)){
    if(is.vector(x)){
      x <- matrix(x, ncol = 1)
    }
    if (nrow(x) == 0) 
      stop("obs has 0 rows")
    if (any(is.na(x))) 
      stop("missing values in obs")
  }
  
  # resp and sumsta recover
  mf <- match.call(expand.dots=FALSE)
  mf <- mf[1]
  mf$formula <- object$formula
  
  mf$data <- training
  
  training <- mf$data
  
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame() )
  mt <- attr(mf, "terms")
  resp <- model.response(mf)
  
  obj <- object$model.rf
  inbag <- matrix(unlist(obj$inbag.counts, use.names=FALSE), ncol=obj$num.trees, byrow=FALSE)
  
  obj[["origNodes"]] <- predict(obj, training, predict.all=TRUE, num.threads=ncores)$predictions
  obj[["origObs"]] <- model.response(mf)
  
  #####################
  
  origObs <- obj$origObs
  origNodes <- obj$origNodes
  
  nodes <- predict(obj, x, predict.all=TRUE, num.threads=ncores )$predictions
  if(is.null(dim(nodes))) nodes <- matrix(nodes, nrow=1)
  ntree <- obj$num.trees
  nobs <- object$model.rf$num.samples
  nnew <- nrow(x)
  weights <- abcrf:::findweights(origNodes, nodes, inbag, 
                                 as.integer(nobs), 
                                 as.integer(nnew), 
                                 as.integer(ntree)) # cpp function call
  
  weights.std <- weights/ntree
  
  post_density <- list()
  
  for(i in seq_len(nnew)) {
    post_density[[i]] <- density(resp, weights = weights.std[, i])
  }
  
  return(post_density)
  
}

sample_density <- function(N, density) {
  
  diffs <- diff(density$x)
  
  # Making sure we have equal increments
  stopifnot(all(abs(diff(density$x) - mean(diff(density$x))) < 1e-9))
  total <- sum(density$y)
  density$y <- density$y / total
  ydistr <- cumsum(density$y)
  yunif <- runif(n)
  indices <- sapply(yunif, function(y) min(which(ydistr > y)))
  x <- density$x[indices]
  
  return(x)
}
