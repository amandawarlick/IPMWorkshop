source("../../R/Rcheck.R")
load.module("glm")
d <- read.jagsdata("epil-data.R")
m <- jags.model("epil2.bug", data=d, n.chains=2)
check.data(m, d)
update(m, 1000)
x <- coda.samples(m, c("alpha0", "alpha.Base", "alpha.Trt", "alpha.BT",
                       "alpha.Age", "alpha.V4", "sigma.b1"), n.iter=5000,
                  thin=5)
source("bench-test1.R")
check.fun()

