source("../../R/Rcheck.R")
d <- read.jagsdata("line-data.R")
m <- jags.model("line.bug", data=d, n.chains=2)
check.data(m, d)
update(m, 1000)
x <- coda.samples(m, c("alpha","beta","sigma","tau"), n.iter=10000)
source("bench-test1.R")
check.fun()
