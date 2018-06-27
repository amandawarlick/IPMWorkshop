#download from here: http://www.openbugs.net/w/Downloads
#follow instructions from here: https://oliviergimenez.github.io/post/run_openbugs_on_mac/

#Resources
#http://www.r-tutor.com/bayesian-statistics/openbugs, http://www.r-tutor.com/content/r-tutorial-ebook
# http://homepage.stat.uiowa.edu/~gwoodwor/BBIText/AppendixBWinbugs.pdf, 
#https://www2.isye.gatech.edu/~brani/isyebayes/bank/shortcoursebugs.pdf, https://www.york.ac.uk/che/pdf/2%20Practical%201%20talk.pdf

library(R2OpenBUGS) 

###from https://oliviergimenez.github.io/post/run_openbugs_on_mac/
data(schools)
nummodel <- function(){
  for (j in 1:J){
    y[j] ~ dnorm (theta[j], tau.y[j])
    theta[j] ~ dnorm (mu.theta, tau.theta)
    tau.y[j] <- pow(sigma.y[j], -2)
    }
  
  mu.theta ~ dnorm (0.0, 1.0E-6)
  tau.theta <- pow(sigma.theta, -2)
  sigma.theta ~ dunif (0, 1000)
}
write.model(nummodel, "nummodel.txt")
model.file1 = paste(getwd(),"nummodel.txt", sep="/")
#file.show("nummodel.txt")

J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list ("J", "y", "sigma.y")

inits <- function(){
  list(theta = rnorm(J, 0, 100),
       mu.theta = rnorm(1, 0, 100),
       sigma.theta = runif(1, 0, 100))
  }

parameters = c("theta", "mu.theta", "sigma.theta")

WINE="/usr/local/Cellar/wine/3.0_1/bin/wine"
WINEPATH="/usr/local/Cellar/wine/3.0_1/bin/winepath"
OpenBUGS.pgm="/Users/amandawarlick/.wine/drive_c/Program\ Files/OpenBUGS/OpenBUGS323/OpenBUGS.exe"

schools.sim <- bugs(data,
                    inits,
                    model.file = model.file1, 
                    parameters=parameters, 
                    n.chains = 3, 
                    n.iter = 100,
                    OpenBUGS.pgm=OpenBUGS.pgm,
                    WINE=WINE,
                    WINEPATH=WINEPATH,
                    useWINE=TRUE,
                    debug = TRUE)



