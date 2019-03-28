rm(list = ls())
library(synlik)

ricker_sl <- synlik(simulator = rickerSimul,
                    summaries = rickerStats,
                    param = c( logR = 3.8, logSigma = log(0.3), logPhi = log(10) ),
                    extraArgs = list("nObs" = 50, "nBurn" = 50)
)

ricker_sl@data <- simulate(ricker_sl, nsim = 1, seed = 54)
ricker_sl@plotFun <- function(input, ...) plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
plot(ricker_sl)

ricker_sl@data <- simulate(ricker_sl, nsim = 1, seed = 54)
ricker_sl@extraArgs$obsData <- ricker_sl@data
tmp <- simulate(ricker_sl, nsim = 2, stats = TRUE)
dim(tmp)
checkNorm(ricker_sl)
slik(ricker_sl, 
     param  = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)),
     nsim   = 1e3)
slice(object = ricker_sl, 
      ranges = list("logR" = seq(3.5, 3.9, by = 0.01),
                    "logPhi" = seq(2, 2.6, by = 0.01),
                    "logSigma" = seq(-2, -0.5, by = 0.02)), 
      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
      nsim = 1000)

slice(object = ricker_sl, 
      ranges = list("logR" = seq(3.5, 3.9, by = 0.01),
                    "logPhi" = seq(2, 2.6, by = 0.01),
                    "logSigma" = seq(-2, -0.5, by = 0.02)), 
      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
      nsim = 1000)

ricker_sl <- smcmc(ricker_sl, 
                   initPar = c(3.2, -1, 2.6),
                   niter = 10, 
                   burn = 3,
                   priorFun = function(input, ...) sum(input), 
                   propCov = diag(c(0.1, 0.1, 0.1))^2, 
                   nsim = 500)
ricker_sl <- continue(ricker_sl, niter = 10)
data(ricker_smcmc)

## If each sim + summary statistics, etc. takes 1 minute
## niter = 20000 nsim = 500 @ each iteration and 3 chains 
20000*500*3 ## 30 million sims = 30 million minutes
3e7/60 ## number of hours
5e5/168 ## number of cores for 1 week
5e5/168/2 ## number of cores for 1 week if each sim takes 30 sec
5e5/168/4 ## number of cores for 1 week if each sim takes 15 sec
5e5/168/6 ## number of cores for 1 week if each sim takes 10 sec


addline1 <- function(parNam, ...) abline(h = ricker_smcmc@param[parNam], lwd = 2, lty = 2, col = 3) 
addline2 <- function(parNam, ...) abline(v = ricker_smcmc@param[parNam], lwd = 2, lty = 2, col = 3)

## Don't need to specify cores, just ncores will try and make a cluster for you--but better to specify and stop it within--need to test on cluster!
system.time({
      slice(object = ricker_sl, 
      ranges = list("logR" = seq(3.5, 3.9, by = 0.02),
                    "logPhi" = seq(2, 2.6, by = 0.02)), 
      pairs = TRUE,
      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
      nsim = 1000, 
      multicore = TRUE, 
      ncores = 2)
})

plot(ricker_smcmc, addPlot1 = "addline1", addPlot2 = "addline2")


check.ncores <- function (sim_secs = 60, niter = 20000, nsim = 500, nchains = 3, target.hrs = 24*7) {
  time_tot = sim_secs*niter*nsim*nchains
  ncores = time_tot/60/60/target.hrs
  print(ncores)
  print(time_tot/60/60)
}

check.ncores(sim_secs = 15, niter = 2000, nsim = 500, nchains = 3, target.hrs = 24*4)

20000*500*3*60 ## 30 million sims = 30 million minutes
1.8e9/60 ## 30 million minutes
3e7/60 ## number of hours
5e5/168 ## number of cores for 1 week
5e5/168/2 ## number of cores for 1 week if each sim takes 30 sec
5e5/168/4 ## number of cores for 1 week if each sim takes 15 sec
5e5/168/6 ## number of cores for 1 week if each sim takes 10 sec
