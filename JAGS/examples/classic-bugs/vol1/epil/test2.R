source("../../R/Rcheck.R")
load.module("glm")
d <- read.jagsdata("epil-data.R")
m <- jags.model("epil3.bug", data=d, n.chains=2)
check.data(m, d)
update(m, 1000)
x <- coda.samples(m, c("alpha0", "alpha.Base", "alpha.Trt", "alpha.BT",
                       "alpha.Age", "alpha.V4", "sigma.b1", "sigma.b"),
                  n.iter=10000, thin=10)
source("bench-test2.R")
check.fun()
