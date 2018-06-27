source("../../R/Rcheck.R")
data <- read.jagsdata("schools-data.R")
inits <- read.jagsdata("schools-inits.R")
load.module("glm")
m <- jags.model("schools.bug", data, inits, n.chains=2)
check.data(m, data)
update(m, 1000)
x <- coda.samples(m, c("beta","gamma","phi","theta"), n.iter=10000, thin=10)
source("bench-test1.R")
check.fun()


