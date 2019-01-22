## Testing binomial probabilities

## Distribution of infectious periods
weeks <- rep(0, 1000)
for (i in 1:1000){
  yes = 0
  while(yes == 0) {
   weeks[i] <- weeks[i] + 1
   yes <- rbinom(1, size = 1, prob = sigma)
  }
}
hist(weeks)
abline(v = mean(weeks), col = "red")

## Distribution of number of incursions on a weekly basis
incursions <- rpois(1000, 1)
hist(incursions)

## Secondary cases generated
incursions <- rnbinom(1000, mu = 1.2, size = 0.4)
hist(incursions)

## Dispersal kernal
shapeM = 0.3484
scaleM = 41.28/100
x<-seq(0,8,0.001)
plot(x,dweibull(x, shape = shapeM, scale = scaleM), ylim=c(0,1.2),
     type="l",bty="l",xlab="Distance to next bite (km)",ylab="Probability Density",
     cex.lab=1.2)
polygon(c(dweibull(x, shape = shapeM, scale = scaleM),rep(0,length(x)))~c(x,rev(x)),col="tomato",border="tomato")
lines(x,dweibull(x, shape = shapeM, scale = scaleM))
