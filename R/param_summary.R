## Testing binomial probabilities
## sigma
get.prob <- function(rate, step) {
  # ## example 1: turn annual waning rate of 0.33 to weekly prob
  # get.prob(rate = 0.33, step = 52)
  # ## example 2: turn annual birth rate of 0.45 to monthly prob
  # get.prob(rate = 0.45, step = 12)
  converted <- (1 + rate)^(1/step) - 1
  return(1 - exp(-converted))
}

sigma <- get.prob(rate = 7/22.3, step = 1)
## Distribution of infectious periods
weeks <- rep(0, 10000)
for (i in 1:10000){
  yes = 0
  while(yes == 0) {
   weeks[i] <- weeks[i] + 1
   yes <- rbinom(1, size = 1, prob = sigma)
  }
}
hist(weeks, border = "white", col = "grey", xlab = "Weeks between bite \n and infectiousness", 
     main = "")
abline(v = mean(weeks), col = "red")

## Distribution of number of incursions on a weekly basis
incursions <- rpois(10000, 1)
hist(incursions)

## Secondary cases generated
secondaries <- rnbinom(10000, mu = 1.1, size = 0.5)
hist(secondaries, include.lowest = TRUE, breaks = seq(0, max(secondaries), by = 1), border = "white", col = "tomato",
     xlab = "Number of secondary cases", 
     main = "")

plot(dnbinom(seq(0, 50, by = 1),  mu = 1.1, size = 0.5), type = "l") ## censored distribution?
     
## Dispersal kernal
shapeM = 0.3484
scaleM = 41.28/100
x <- seq(0, 8, 0.001)
plot(x, dweibull(x, shape = shapeM, scale = scaleM), ylim=c(0,1.2),
     type = "l", bty = "l", xlab = "Distance to next bite (km)", ylab = "Probability Density",
     cex.lab = 1.2)
polygon(c(dweibull(x, shape = shapeM, scale = scaleM),rep(0,length(x)))~c(x,rev(x)),col="tomato",border="tomato")
lines(x,dweibull(x, shape = shapeM, scale = scaleM))
