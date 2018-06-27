source("../../R/Rcheck.R")
data <- read.jagsdata("jaw-data.R")
data$PI <- NULL
inits <- read.jagsdata("jaw-inits.R")
inits$data2 <- NULL #Not used in this model
m <- jags.model("jaw-linear.bug", data, inits, n.chains=2)
check.data(m, data)
update(m, 1000)
load.module("dic")
x <- coda.samples(m,
                  c("beta0.uncentred","beta1","Sigma2","mu","RSS","deviance"),
                  n.iter=10000)
source("bench-test2.R")
check.fun()

