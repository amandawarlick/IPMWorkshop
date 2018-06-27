source("../../R/Rcheck.R")
data <- read.jagsdata("schools-data.R")
inits <- read.jagsdata("schools-inits.R")
load.module("glm")
##set.factory("glm::REScaledWishart", FALSE, type="sampler")
set.factory("glm::REScaledWishart", TRUE, type="sampler")
m <- jags.model("schools2.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("beta","gamma","phi","theta"), n.iter=10000, thin=10)
source("bench-test2.R")
check.fun()


