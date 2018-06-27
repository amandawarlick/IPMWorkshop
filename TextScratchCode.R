#page 19
reps <- 10^6
sample.means <- rep(NA, reps)

mu <- 65
sigma <- 5

for(i in 1:reps){
  sample.means[i] <- mean(rnorm(n = 10, mean = mu, sd = sigma))
}

x <- sample.means

#plot the distribution of one million samples of 10
qplot(sample.means)
par(mfrow = c(1, 2), las = 1)
hist(x, col = "grey", main = "", xlab = "Body length (cm)", las = 1)
abline(v = mu, lwd = 3, col = "red")
abline(v = mean(x), lwd = 3, col = "blue")

#Chapter 3

##Peregrin falcons; linear predictor is cubic binomial function

data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014) {
  #generating values of time covariate
  year <- 1:n 
  #signal part of GLM
  log.expected.count <- alpha + beta1*year + beta2*year^2 + beta3*year^3 #Q:why is this additive and not just yr^3?
  expected.count <- exp(log.expected.count)
  #poisson noise around expected counts
  C <- rpois(n = n, lambda = expected.count) 
  #plot simulated data; Q: dots are supposed to be "observed data", which is the poisson noise C
  plot(year, C, type = "b", col = "black", las = 1)
  lines(year, expected.count, col = "red")
  return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3,
              year = year, expected.count = expected.count, C = C))
}

data <- data.fn()

#plot(data$year, data$C)

#Analyze this data set/function using R

fm <- glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data)
summary(fm)

#Set-up and analyze in BUGS
library(R2OpenBUGS)
WINE="/usr/local/Cellar/wine/3.0_1/bin/wine"
WINEPATH="/usr/local/Cellar/wine/3.0_1/bin/winepath"
OpenBUGS.pgm="/Users/amandawarlick/.wine/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

#specify model in BUGS
sink("GLM_Poisson.txt")
cat("
    model {
#priors
    alpha ~ dunif(-20, 20)
    beta1 ~ dunif(-10, 10)
    beta2 ~ dunif(-10, 10)
    beta3 ~ dunif(-10, 10)
#likelihood
  for(i in 1:n {
C[i]} ~ dpois(lambda[i]) #distribution of random element
log(lambda[i]) <- log.lambda[i] #link function
log(lambda[i] <- alpha + beta1*year[i] + 
        beta2*pow(year[i], 2) + beta3*pow(year[i], 3) #linear predictor
    }
}
", fill = TRUE) 
sink()

#bundle up the data
#win.data <- list(C = data$C, n = length(data$C), year = data$year)

#initial values - must define 1
inits <- function () list(alpha = runif(1, -2, 2),
                          beta1 = runif(1, -3, 3)) 

#parameters
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

#MCMC
ni <- 2000
nt <- 2
nb <- 1000
nc <- 3

#center/standardize covariates (year^3 is too large otherwise)
mean.year <- mean(data$year)
sd.year <- sd(data$year)
win.data <- list(C = data$C, n = length(data$C), year = (data$year - mean.year)/sd.year)

out <- bugs(data = win.data, inits = inits, 
            parameters.to.save = params,
            model.file = "GLM_Poisson.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni,
            n.burnin = nb, OpenBUGS.pgm=OpenBUGS.pgm,
            WINE=WINE,
            WINEPATH=WINEPATH,
            useWINE=TRUE,
            debug = TRUE)

#Chapter 12: N-mix models

R <- 200
T <- 3

y <- array(dim = c(R, T))
N <- rpois(n = R, lambda = 2)

