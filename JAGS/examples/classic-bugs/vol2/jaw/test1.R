source("../../R/Rcheck.R")
data <- read.jagsdata("jaw-data.R")
data$PI <- NULL
inits <- read.jagsdata("jaw-inits.R")
inits[c("beta1","beta2")] <- NULL #Not used in this model
m <- jags.model("jaw-constant.bug", data, inits, n.chains=2)
check.data(m, data)
update(m, 1000)
load.module("dic")
x <- coda.samples(m, c("beta0", "mu", "Sigma2", "RSS", "deviance"),
                  n.iter=10000)
source("bench-test1.R")
check.fun()

