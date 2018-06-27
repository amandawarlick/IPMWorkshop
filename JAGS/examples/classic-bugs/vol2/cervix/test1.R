load.module("glm")
data <- read.jagsdata("cervix-data.R")
m <- jags.model("cervix.bug", data=data,
                inits=read.jagsdata("cervix-inits.R"), n.chain=2)
check.data(m, data)
update(m, 1000)
x <- coda.samples(m, c("beta0C","beta","phi","gamma1","gamma2"),
                  n.iter=10000, thin=10)
source("bench-test1.R")
check.fun()


