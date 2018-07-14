########################################################################
#
# State-space models 
#
###########################################################################
#
# IPM workshop Aberdeen June 2018
#
###########################################################################
#
# Based on the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc K?ry & Michael Schaub (2012, Academic Press)
#
###########################################################################
#
# MS, 12 September 2017
#
###########################################################################


# Load packages
library(jagsUI)


###########################################################################
#
# Example 1: state-space models for a stochastic Gaussian linear system
#
###########################################################################

# Data simulation for a random walk, starting from 0

ss.sim <- function(sigma.proc = 5, sigma.obs = 10, T = 100, n1 = 0){
  true <- observed <- numeric()
  true[1] <- n1
  for (i in 1:(T-1)){
    true[i+1] <- rnorm(1, true[i], sigma.proc)
  }
  for (i in 1:T){
    observed[i] <- rnorm(1, true[i], sigma.obs)
  }
  
  plot(true, type = "l", ylim = range(c(true, observed)), xlab = "Time", ylab = "x", las = 1, col = "blue")
  lines(observed, col = "red")
  legend("topright", legend = c("True", "Observed"), col = c("blue", "red"), lty = c(1,1), bty = "n")
  return(list(observed = observed, true = true))
}

z <- ss.sim(sigma.proc = 7, sigma.obs = 10, T = 50)


# Specify model in BUGS language
cat(file = "ssm.jags", "
    model { 
    
    # Priors and constraints
    x[1] ~ dnorm(0, 0.01)            # Prior for initial population size
    sigma1 ~ dunif(0, 100)           # Prior for sd of state process
    tau1 <- pow(sigma1, -2)
    
    sigma2 ~ dunif(0, 500)           # Prior for sd of observation process; no specific reason it is not 100 as above
    tau2 <- pow(sigma2, -2)
    
    # Likelihood
    # State process: loop over values 2-50 bc Markovian process, T = 1 = prior?
    for (t in 1:(T-1)){
    x[t+1] ~ dnorm(x[t], tau1) 
    }
    # Observation process: looping over time
    for (t in 1:T) {
    y[t] ~ dnorm(x[t], tau2)        #In BUGS/JAGS, dnorm is mean/precision, not variance
    }
    }
    ")

# Bundle data
jags.data <- list(y = z$observed, T = length(z$observed))

# Initial values
inits <- function(){list(sigma1 = runif(1, 0, 5))} 

# Parameters monitored
parameters <- c("x", "sigma1", "sigma2")

# MCMC settings
ni <- 5000
nt <- 3
nb <- 1000
nc <- 3

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)



# Define function to draw a graph to summarize results
graph.ssm <- function(ssm, x, y){
  fitted <- lower <- upper <- numeric()
  n.years <- length(y)
  for (i in 1:n.years){
    fitted[i] <- mean(ssm$sims.list$x[,i])
    lower[i] <- quantile(ssm$sims.list$x[,i], 0.025)
    upper[i] <- quantile(ssm$sims.list$x[,i], 0.975)}
  m1 <- min(c(y, fitted, x, lower))
  m2 <- max(c(y, fitted, x, upper))
  par(mar = c(4.5, 4, 1, 1), cex = 1.2)
  plot(0, 0, ylim = c(m1, m2), xlim = c(0.5, n.years), ylab = "X", xlab = "Time", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
  axis(2, las = 1)
  axis(1)
  polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  points(x, type = "l", col = "blue", lwd = 2)
  points(y, type = "l", col = "red", lwd = 2)
  points(fitted, type = "l", col = "black", lwd = 2)
  legend("topright", legend = c("Truth", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("blue", "red", "black"), bty = "n", cex = 1)
}

# Execute function: Produce figure 
graph.ssm(ssm, z$true, z$observed)

hist(ssm$sims.list$sigma1, nclass = 20)



###########################################################################
#
# Example 2: state-space models for movement data
#
###########################################################################

# Data simulation function

ss.sim <- function(sigmaX.proc = 5, sigmaY.proc = 5, sigma.obs = 10, T = 100, x1 = 0, y1 = 0){
  coord.true <- coord.obs <- matrix(NA, ncol = 2, nrow = T)
  coord.true[1,1] <- x1
  coord.true[1,2] <- y1
  
  for (i in 1:(T-1)){
    coord.true[i+1,1] <- rnorm(1, coord.true[i,1], sigmaX.proc)
    coord.true[i+1,2] <- rnorm(1, coord.true[i,2], sigmaY.proc)
  }
  for (i in 1:T){
    coord.obs[i,1] <- rnorm(1, coord.true[i,1], sigma.obs)
    coord.obs[i,2] <- rnorm(1, coord.true[i,2], sigma.obs)
  }
  
  plot(x = coord.true[,1], y = coord.true[,2], pch = 16, xlim = range(c(coord.true[,1], coord.obs[,1])), ylim = range(c(coord.true[,2], coord.obs[,2])), xlab = "X", ylab = "Y", col = "blue")
  for (i in 1:(T-1)){
    segments(coord.true[i+1,1], coord.true[i+1,2], coord.true[i,1], coord.true[i,2], col = "blue")
  }
  points(x = coord.obs[,1], y = coord.obs[,2], pch = 16, , col = "red")
  for (i in 1:(T-1)){
    segments(coord.obs[i+1,1], coord.obs[i+1,2], coord.obs[i,1], coord.obs[i,2], col = "red")
  }
  points(x = coord.true[1,1], y = coord.true[1,2], pch = 15, cex = 1.5, col = "blue")
  points(x = coord.obs[1,1], y = coord.obs[1,2], pch = 15, cex = 1.5, col = "red")
  legend("topright", legend = c("Truth", "Observed"), lty = c(1, 1), lwd = c(2, 2), col = c("blue", "red"), bty = "n", cex = 1)
  return(list(coord.true = coord.true, coord.obs = coord.obs))
}

z <- ss.sim(sigmaX.proc = 10, sigmaY.proc = 12, sigma.obs = 5, T = 30)


# Specify model in BUGS language
cat(file = "ssm.jags", "
    model { 
    
    # Priors and constraints
    x[1] ~ dnorm(0, 0.01)            # Prior for initial population size
    sigmaX ~ dunif(0, 100)           # Prior for sd of state process
    tauX <- pow(sigmaX, -2)
    
    y[1] ~ dnorm(0, 0.01)            # Prior for initial population size
    sigmaY ~ dunif(0, 100)           # Prior for sd of state process
    tauY <- pow(sigmaY, -2)
    
    sigma.obs ~ dunif(0, 100)           # Prior for sd of observation process
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    x[t+1] ~ dnorm(x[t], tauX)
    y[t+1] ~ dnorm(y[t], tauY)   
    }
    # Observation process
    for (t in 1:T) {
    x.obs[t] ~ dnorm(x[t], tau.obs)
    y.obs[t] ~ dnorm(y[t], tau.obs)
    }
    }
    ")



# Bundle data
jags.data <- list(x.obs = z$coord.obs[,1], y.obs = z$coord.obs[,2], T = nrow(z$coord.obs))

# Initial values
inits <- function(){list(sigmaX = runif(1, 0, 20))} 

# Parameters monitored
parameters <- c("x", "y", "sigmaX", "sigmaY", "sigma.obs")

# MCMC settings
ni <- 5000
nt <- 6
nb <- 2500
nc <- 3

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)



# Define function to draw a graph to summarize results
graph.ssm <- function(ssm, coord.true, coord.obs, err = NA){
  T <- nrow(coord.true)
  fitted <- lower <- upper <- matrix(NA, ncol = 2, nrow = T)
  for (i in 1:T){
    fitted[i,1] <- mean(ssm$sims.list$x[,i])
    fitted[i,2] <- mean(ssm$sims.list$y[,i])
  }
  x.m <- range(c(coord.true[,1], fitted[,1], coord.obs[,1], ssm$sims.list$x))
  y.m <- range(c(coord.true[,2], fitted[,2], coord.obs[,2], ssm$sims.list$y))
  par(mar = c(4.5, 4, 1, 1), cex = 1.2)
  
  plot(x = coord.true[,1], y = coord.true[,2], pch = 16, xlim = x.m, ylim = y.m, xlab = "X", ylab = "Y")
  if (is.na(err)==FALSE){
    u <- nrow(ssm$sims.list$x)
    for (k in 1:u){
      for (i in 1:(T-1)){
        segments(ssm$sims.list$x[k,i+1], ssm$sims.list$y[k,i+1], ssm$sims.list$x[k,i], ssm$sims.list$y[k,i], col = "lightgrey", lwd = 0.5)
      }
    }
  }
  
  points(x = coord.true[,1], y = coord.true[,2], pch = 16, col = "blue")
  for (i in 1:(T-1)){
    segments(coord.true[i+1,1], coord.true[i+1,2], coord.true[i,1], coord.true[i,2], col = "blue")
  }
  points(x = coord.obs[,1], y = coord.obs[,2], pch = 16, , col = "red")
  for (i in 1:(T-1)){
    segments(coord.obs[i+1,1], coord.obs[i+1,2], coord.obs[i,1], coord.obs[i,2], col = "red")
  }
  points(x = fitted[,1], y = fitted[,2], pch = 16, col = "black")
  for (i in 1:(T-1)){
    segments(fitted[i+1,1], fitted[i+1,2], fitted[i,1], fitted[i,2], col = "black")
  }
  points(x = coord.true[1,1], y = coord.true[1,2], pch = 15, cex = 1.5, col = "blue")
  points(x = coord.obs[1,1], y = coord.obs[1,2], pch = 15, cex = 1.5, col = "red")
  points(x = fitted[1,1], y = fitted[1,2], pch = 15, cex = 1.5, col = "black")
  
  legend("topright", legend = c("Truth", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("blue", "red", "black"), bty = "n", cex = 1)
}

# Execute function: produce figure 
graph.ssm(ssm, z$coord.true, z$coord.obs)                  # without measure of uncertainty
graph.ssm(ssm, z$coord.true, z$coord.obs, err = "y")        # with measure of uncertainty





###########################################################################
#
# State-space models applied for population count data (chapter 5 in the BPA book)
#
###########################################################################


# 5.2. A simple model
n.years <- 25           # Number of years
N1 <- 30                # Initial population size
mean.lambda <- 1.02     # Mean annual population growth rate
sigma2.lambda <- 0.02   # Process (temporal) variation of the growth rate
sigma2.y <- 20          # Variance of the observation error

y <- N <- numeric(n.years)
N[1] <- N1
lambda <- rnorm(n.years-1, mean.lambda, sqrt(sigma2.lambda))
for (t in 1:(n.years-1)){
  N[t+1] <- N[t] * lambda[t]
}

for (t in 1:n.years){
  y[t] <- rnorm(1, N[t], sqrt(sigma2.y))
}


# Specify model in BUGS language
cat(file = "ssm.jags", "
    model { 

    # Priors and constraints; unif bc has to be positive (not dnorm)
    N.est[1] ~ dunif(0, 500)            # Prior for initial population size
    mean.lambda ~ dunif(0, 10)          # Prior for mean growth rate
    sigma.proc ~ dunif(0, 10)           # Prior for SD of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 100)           # Prior for SD of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    ## State process
    for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mean.lambda, tau.proc) 
    N.est[t+1] <- N.est[t] * lambda[t] 
    }
    ## Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(N.est[t], tau.obs)
    }
    }
    ")

# Bundle data
jags.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 20, 40), rep(NA, (n.years-1))))} 

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


# Define function to draw a graph to summarize results
graph.ssm <- function(ssm, N, y){
  fitted <- lower <- upper <- numeric()
  n.years <- length(y)
  for (i in 1:n.years){
    fitted[i] <- mean(ssm$sims.list$N.est[,i])
    lower[i] <- quantile(ssm$sims.list$N.est[,i], 0.025)
    upper[i] <- quantile(ssm$sims.list$N.est[,i], 0.975)}
  m1 <- min(c(y, fitted, N, lower))
  m2 <- max(c(y, fitted, N, upper))
  par(mar = c(4.5, 4, 1, 1), cex = 1.2)
  plot(0, 0, ylim = c(m1, m2), xlim = c(0.5, n.years), ylab = "Population size", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
  axis(2, las = 1)
  axis(1, at = seq(0, n.years, 5), labels = seq(0, n.years, 5))
  axis(1, at = 0:n.years, labels = rep("", n.years + 1), tcl = -0.25)
  polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  points(N, type = "l", col = "red", lwd = 2)
  points(y, type = "l", col = "black", lwd = 2)
  points(fitted, type = "l", col = "blue", lwd = 2)
  legend(x = 1, y = m2, legend = c("True", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"), bty = "n", cex = 1)
}

# Execute function: Produce figure 
graph.ssm(ssm, N, y)



###########################################################################
#
# 5.3. Systematic bias in the observation process
#
###########################################################################

n.years <- 25  # Number of years
N <- rep(50, n.years) 

p <- 0.7
y <- numeric(n.years)
for (t in 1:n.years){
  y[t] <- rbinom(1, N[t], p)
}
y

# Bundle data
jags.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 30, 60), rep(NA, (n.years-1))))}

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Summarize posteriors
print(ssm, digits = 3)

# Produce figure
graph.ssm(ssm, N, y)



# Trend in detection probability

n.years <- 25  # Number of years
N <- rep(50, n.years)

lp <- -0.5 + 0.1*(1:n.years)  # Increasing trend of logit p
p <- plogis(lp)
y <- numeric(n.years)
for (t in 1:n.years){
  y[t] <- rbinom(1, N[t], p[t])
}

# Bundle data
jags.data <- list(y = y, T = n.years)

# Call JAGS from R (BRT <1 min)
ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# Produce figure
graph.ssm(ssm, N, y)
points(N*p, col = "black", type = "l", lwd = 2, lty = 2)
legend(x = 1, y = 45.5, legend = "Np", lwd = 2, col = "black", lty = 2, bty = "n")


###########################################################################
#
# 5.4. Real example: House martin population counts in the village of Magden
#
###########################################################################


# Specify model in BUGS language
cat(file = "ssm.jags", "
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
    }
    }
    ")

# House martin population data from Magden
pyears <- 7 # Number of future years with predictions
hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears))
year <- 1990:(2009 + pyears)

# Bundle data
jags.data <- list(y = log(hm), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), 
                         logN.est = c(rnorm(1, 5.6, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 3 min)
hm.ssm <- jags(jags.data, inits, parameters, "ssm.jags", n.chains = nc, n.thin = nt, 
               n.iter = ni, n.burnin = nb, parallel = T)

# Summarize posteriors
print(hm.ssm, digits = 3)


# Draw figure
fitted <- lower <- upper <- numeric()
year <- 1990:2016
n.years <- length(hm)
for (i in 1:n.years){
  fitted[i] <- mean(hm.ssm$sims.list$N.est[,i])
  lower[i] <- quantile(hm.ssm$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(hm.ssm$sims.list$N.est[,i], 0.975)}
m1 <- min(c(fitted, hm, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)

# Probability of N(2016) < N(2009)
mean(hm.ssm$sims.list$N.est[,27] < hm.ssm$mean$N.est[20])

# How well are our predictions? Add data from the years 2010 to 2014
points(y = c(180, 202, 224, 203, 194), x = c(21, 22, 23, 24, 25), col = "black", pch = 16)


###########################################################################
#
# Exercises
#
###########################################################################


# Exercise 1: From year 1 to 10 a different data sampling protocol than in years 11 to 20 has ben used. 
#Adapt the model to account for possible different observation errors

#solution 1a: break up observation process loop into two different time periods or
#solution 1b: use categorical covariate

cat(file = "ssm.jags", "
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)

    sigma.obs1 ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs1 <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)

    sigma.obs2 ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.ob2 <- pow(sigma.obs1, 2)
    tau.obs <- pow(sigma.obs, -2)
  

    # Likelihood

    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs[survey[t]]) #add covariate for survey type
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
    }
    }
    ")

#############################

# Exercise 2: model the population growth rate of the house martins as a function of rainfall during the summer.

# Data on precipitation 1990 - 2009
May <- c(33.1, 42.8, 33.7, 89.8, 172.7, 180.7, 89.2, 59.1, 49.5, 141.3, 54.6, 117.8, 142.5, 80.4, 50.9, 98.3, 140.9, 160.2, 74.6, 86.7)
June <- c(155.0, 138.8, 124.4, 89.4, 49.7, 45.5, 77.9, 133.2, 75.0, 126.4, 67.4, 120.4, 79.3, 21.0, 70.1, 84.2, 76.3, 183.1, 34.3, 65.9)
rain <- May + June

cat(file = "ssm.jags", "
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
    }
    }
    ")


#############################

# Exercise 3: Assume that we know (or have a predction) precipitation in the future: 
#what is our prediction of the house martin population?

# Data from 1990 to 2016
May <- c(33.1, 42.8, 33.7, 89.8, 172.7, 180.7, 89.2, 59.1, 49.5, 141.3, 54.6, 117.8, 142.5, 80.4, 50.9, 98.3, 140.9, 160.2, 74.6, 86.7, 140.5, 48.4, 66.5, 99.4, 70, 92.1, 154.7)
June <- c(155.0, 138.8, 124.4, 89.4, 49.7, 45.5, 77.9, 133.2, 75.0, 126.4, 67.4, 120.4, 79.3, 21.0, 70.1, 84.2, 76.3, 183.1, 34.3, 65.9, 65.4, 96.9, 134.3, 64.5, 59.4, 94.5, 176.2)
rain <- May + June

cat(file = "ssm.jags", "
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
    }
    }
    ")


#############################

# Exercise 4: Assume that we have an estimate (mean and measure of uncertainty) for the precipitation in the future: 
#what is our prediction of the house martin population? 

# Data from 1990 to 2016
May <- c(33.1, 42.8, 33.7, 89.8, 172.7, 180.7, 89.2, 59.1, 49.5, 141.3, 54.6, 117.8, 142.5, 80.4, 50.9, 98.3, 140.9, 160.2, 74.6, 86.7, 140.5, 48.4, 66.5, 99.4, 70, 92.1, 154.7)
June <- c(155.0, 138.8, 124.4, 89.4, 49.7, 45.5, 77.9, 133.2, 75.0, 126.4, 67.4, 120.4, 79.3, 21.0, 70.1, 84.2, 76.3, 183.1, 34.3, 65.9, 65.4, 96.9, 134.3, 64.5, 59.4, 94.5, 176.2)
rain <- May + June
# Accuracy of the predictions in the future
sd.rain <- c(5, 5, 5, 10, 10, 20, 20)

cat(file = "ssm.jags", "
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(0, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
    }
    }
    ")
