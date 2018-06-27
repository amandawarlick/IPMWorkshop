source("../../R/Rcheck.R")
d <- read.jagsdata("equiv-data.R")
d$T <- NULL #Variable T not used in this model
m <- jags.model("equiv.bug", d, n.chains=2)
check.data(m, d)
update(m, 1000)
x <- coda.samples(m, c("theta","equivalence","sigma"), n.iter=10000)
source("bench-test1.R")
check.fun()
