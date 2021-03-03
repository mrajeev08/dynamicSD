# Constrain by max size because this is what slows it down!

sim.bp <- function(R0 = 1.1, k = 1, nsim = 5, max_size = 1000) {
  progen <- chain_size <- 1
  length <- 0
  sizes <- lengths <- vector(length = nsim)
  
  for(i in 1:nsim){
    counter <- 0
    while (progen > 0 & chain_size < max_size){
      vec <- rnbinom(n = progen, size = k, mu = R0)
      progen <- sum(vec)
      chain_size <- chain_size + progen
      counter <- counter + 1
      chain_length <- ifelse(progen > 0, counter, counter -1)
    }
    sizes[i] <- chain_size
    lengths[i] <- chain_length
    progen <- chain_size <- 1 # set for next loop
  }
  return(list(sizes, lengths))
}

check <- sim.bp(nsim = 1)

# profile with chain_sim as well!