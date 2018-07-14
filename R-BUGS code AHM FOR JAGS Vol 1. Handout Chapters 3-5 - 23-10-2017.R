
################################################################
##  Code for part of the chapters 3 -- 5 in the book
##
##  Applied hierarchical modeling in ecology
##  Modeling distribution, abundance and species richness using R and BUGS
##  Volume 1: Prelude and Static models, 
##  Marc Kéry & J. Andy Royle
##  Academic Press, 2016
##
##
##  This document contains part of Chapters 3–5
##
##  3 September 2017, 23 October 2017
##
################################################################


####################################################
#   3. Linear models, generalised linear models (GLMs) 
#   and random effects models: 
#   the components of hierarchical models
####################################################


### 3.1 Introduction

# Define data
pop <- factor(c(rep("Navarra", 3), rep("Aragon", 3), rep("Catalonia", 3)), levels = c("Navarra", "Aragon", "Catalonia"))         # Population
wing <- c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9) # Wing span
body <- c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2) # Body length
sex <- factor(c("M","F","M","F","M","F","M","F","M"), levels = c("M", "F"))
mites <- c(0, 3, 2, 1, 0, 7, 0, 9, 6)      # Number of ectoparasites
color <- c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57) # Color intensity
damage <- c(0,2,0,0,4,2,1,0,1)                 # Number of wings damaged

cbind(pop, sex, wing, body, mites, color, damage) # Look at data

str(pop)

par(mfrow = c(1, 3), cex = 1.2)
colorM <- c("red", "red", "blue", "green", "green")  # Pop color code males
colorF <- c("red", "blue", "blue", "green", "green") # Pop color code females
plot(body[sex == "M"], wing[sex == "M"], col =colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Wing span")
points(body[sex == "F"], wing[sex == "F"], col =colorF, pch = 16)
text(6.8, 13.8, "A", cex = 1.5)
plot(body[sex == "M"], mites[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(0, 10), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Parasite load")
points(body[sex == "F"], mites[sex == "F"], col = colorF, pch = 16)
text(6.8, 9.7, "B", cex = 1.5)
plot(body[sex == "M"], damage[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(0, 4), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Damaged wings")
points(body[sex == "F"], damage[sex == "F"], col = colorF, pch = 16)
text(6.8, 3.9, "C", cex = 1.5)


### 3.2 Linear models

      
### 3.2.1 Linear models with main effects of one factor and one continuous covariate

summary(fm1 <- lm(wing ~ pop + body))

summary(fm2 <- lm(wing ~ pop-1 + body))

cbind(model.matrix(~pop+body) %*% fm1$coef, predict(fm1)) # Compare two solutions

model.matrix(~ pop + body) # Effects parameterisation

model.matrix(~ pop-1 + body) # Means parameterization 

par(mfrow = c(1, 3), mar = c(5,4,2,2), cex = 1.2, cex.main = 1)
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Wing span")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm2)[1], coef(fm2)[4], col = "red", lwd = 2)
abline(coef(fm2)[2], coef(fm2)[4], col = "blue", lwd = 2)
abline(coef(fm2)[3], coef(fm2)[4], col = "green", lwd = 2)
text(6.8, 14, "A", cex = 1.5)


### 3.2.2 Linear models with interaction between one factor and one continuous covariate

model.matrix(~ pop*body)  # Effects parameterisation

model.matrix(~ pop*body-1-body)  # Means parameterisation
# Output slightly trimmed

summary(fm3 <- lm(wing ~ pop*body-1-body))
Call:
lm(formula = wing ~ pop * body - 1 - body)

# Plot
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm3)[1], coef(fm3)[4], col = "red", lwd = 2)
abline(coef(fm3)[2], coef(fm3)[5], col = "blue", lwd = 2)
abline(coef(fm3)[3], coef(fm3)[6], col = "green", lwd = 2)
text(6.8, 14, "B", cex = 1.5)


# Create new design matrix
(DM0 <- model.matrix(~ pop*body-1-body)) # Original DM for means param
DM0[7:9,5] <- DM0[7:9,6]                 # Combine slopes for Ar and Cat
(DM1 <- DM0[, -6])                       # Delete former slope column for Cat

# Fit model with partial interaction
summary(fm4 <- lm(wing ~ DM1-1))


# Do significance test
anova(fm3, fm4)             # F test between two models

# Plot
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm4)[1], coef(fm4)[4], col = "red", lwd = 2)
abline(coef(fm4)[2], coef(fm4)[5], col = "blue", lwd = 2)
abline(coef(fm4)[3], coef(fm4)[5], col = "green", lwd = 2)
text(6.8, 14, "C", cex = 1.5)



### 3.3 Generalised linear models (GLMs)

### 3.3.1 Poisson generalised linear model (GLM) for unbounded counts

summary(fm10 <- glm(mites ~ pop-1 + body, family = poisson))


### 3.3.5 Bernoulli GLM: logistic regression for a binary response

presence <- ifelse(mites > 0, 1, 0)  # convert abundance to presence/absence
summary(fm11 <- glm(presence ~ pop-1 + body, family = binomial))


### 3.3.7 Binomial GLM: logistic regression for bounded counts

summary(fm12 <- glm(cbind(damage, 4-damage) ~ pop + body -1, family = binomial))


### 3.4 Random effects (mixed) models

### 3.4.1 Random effects for a normal data distribution: normal-normal generalised linear mixed model (GLMM)

# Plot data without distinguishing sex
plot(body, wing, col = rep(c("red", "blue", "green"), each = 3), xlim = c(6.5, 9.5), ylim = c(10, 14), cex = 1.5, lwd = 2, frame.plot = FALSE, las = 1, pch = 16, xlab = "Body length", ylab = "Wing span")

summary(lm <- lm(wing ~ pop-1 + body))     # Same as fm2

library(lme4)
summary(lmm1 <- lmer(wing ~ (1|pop) + body))  # Fit the model
ranef(lmm1)                                   # Print random effects


alpha_j <- fixef(lmm1)[1]+ranef(lmm1)$pop[,1]
cbind(fixed = coef(lm)[1:3], random = alpha_j)

par(lwd = 3)
abline(lm$coef[1], lm$coef[4], col = "red", lty = 2)
abline(lm$coef[2], lm$coef[4], col = "blue", lty = 2)
abline(lm$coef[3], lm$coef[4], col = "green", lty = 2)
abline(alpha_j[1], fixef(lmm1)[2], col = "red")
abline(alpha_j[2], fixef(lmm1)[2], col = "blue")
abline(alpha_j[3], fixef(lmm1)[2], col = "green")
abline(fixef(lmm1), col = "black")
legend(6.5, 14, c("Catalonia", "Aragon", "Navarra"), col=c("blue", "green", "red"), lty = 1, pch = 16, bty = "n", cex = 1.5)

summary(lmm2 <- lmer(wing ~ body + (1|pop) + (0+body|pop)))


### 3.4.2 Random effects for a Poisson data distribution: normal-Poisson generalised linear mixed model (GLMM)

summary(glmm <- glmer(mites ~ body + (1|pop), family = poisson))



#######################################
#  4. Introduction to data simulation
#######################################

### 4.1 What do we mean by data simulation and why is it so tremendously useful ?

### 4.2 Generation of a typical point count data set

### 4.3 Packaging everything in a function
# Function definition with set of default values
data.fn <- function(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, 
                    mean.detection = 0.3, alpha1 = 1, alpha2 = -3, alpha3 = 0, show.plot = TRUE){
#
# Function to simulate point counts replicated at M sites during J occasions.
# Population closure is assumed for each site.
# Expected abundance may be affected by elevation (elev), 
# forest cover (forest) and their interaction.
# Expected detection probability may be affected by elevation, 
# wind speed (wind) and their interaction.
# Function arguments:
#     M: Number of spatial replicates (sites)
#     J: Number of temporal replicates (occasions)
#     mean.lambda: Mean abundance at value 0 of abundance covariates
#     beta1: Main effect of elevation on abundance
#     beta2: Main effect of forest cover on abundance
#     beta3: Interaction effect on abundance of elevation and forest cover
#     mean.detection: Mean detection prob. at value 0 of detection covariates
#     alpha1: Main effect of elevation on detection probability
#     alpha2: Main effect of wind speed on detection probability
#     alpha3: Interaction effect on detection of elevation and wind speed
#     show.plot: if TRUE, plots of the data will be displayed; 
#        set to FALSE if you are running simulations.

# Create covariates
elev <- runif(n = M, -1, 1)                         # Scaled elevation
forest <- runif(n = M, -1, 1)                       # Scaled forest cover
wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed

# Model for abundance
beta0 <- log(mean.lambda)               # Mean abundance on link scale
lambda <- exp(beta0 + beta1*elev + beta2*forest + beta3*elev*forest)
N <- rpois(n = M, lambda = lambda)      # Realised abundance
Ntotal <- sum(N)                        # Total abundance (all sites)
psi.true <- mean(N>0)                   # True occupancy in sample

# Plots
if(show.plot){
par(mfrow = c(2, 2), cex.main = 1)
devAskNewPage(ask = TRUE)
curve(exp(beta0 + beta1*x), -1, 1, col = "red", main = "Relationship lambda-elevation \nat average forest cover", frame.plot = F, xlab = "Scaled elevation")
plot(elev, lambda, xlab = "Scaled elevation", main = "Relationship lambda-elevation \nat observed forest cover", frame.plot = F)
curve(exp(beta0 + beta2*x), -1, 1, col = "red", main = "Relationship lambda-forest \ncover at average elevation", xlab = "Scaled forest cover", frame.plot = F)
plot(forest, lambda, xlab = "Scaled forest cover", main = "Relationship lambda-forest cover \nat observed elevation", frame.plot = F)
}

# Model for observations
alpha0 <- qlogis(mean.detection)        # mean detection on link scale
p <- plogis(alpha0 + alpha1*elev + alpha2*wind + alpha3*elev*wind)
C_jags <- matrix(NA, nrow = M, ncol = J)     # Prepare matrix for counts
for (i in 1:J){                         # Generate counts by survey
   C[,i] <- rbinom(n = M, size = N, prob = p[,i])
}
summaxC <- sum(apply(C,1,max))          # Sum of max counts (all sites)
psi.obs <- mean(apply(C,1,max)>0)       # Observed occupancy in sample

# More plots
if(show.plot){
par(mfrow = c(2, 2))
curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", main = "Relationship p-elevation \nat average wind speed", xlab = "Scaled elevation", frame.plot = F)
matplot(elev, p, xlab = "Scaled elevation", main = "Relationship p-elevation\n at observed wind speed", pch = "*", frame.plot = F)
curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", main = "Relationship p-wind speed \n at average elevation", xlab = "Scaled wind speed", frame.plot = F)
matplot(wind, p, xlab = "Scaled wind speed", main = "Relationship p-wind speed \nat observed elevation", pch = "*", frame.plot = F)

matplot(elev, C, xlab = "Scaled elevation", main = "Relationship counts and elevation", pch = "*", frame.plot = F)
matplot(forest, C, xlab = "Scaled forest cover", main = "Relationship counts and forest cover", pch = "*", frame.plot = F)
matplot(wind, C, xlab = "Scaled wind speed", main = "Relationship counts and wind speed", pch = "*", frame.plot = F)
desc <- paste('Counts at', M, 'sites during', J, 'surveys')
hist(C, main = desc, breaks = 50, col = "grey")
}

# Output
return(list(M = M, J = J, mean.lambda = mean.lambda, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3, mean.detection = mean.detection, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, elev = elev, forest = forest, wind = wind, lambda = lambda, N = N, p = p, C = C, Ntotal = Ntotal, psi.true = psi.true, summaxC = summaxC, psi.obs = psi.obs))
}


data.fn()                  # Execute function with default arguments

data.fn(show.plot = FALSE) # same, without plots
data.fn(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 0.3, alpha1 = 1, alpha2 = -3, alpha3 = 0) # Explicit defaults
data <- data.fn()          # Assign results to an object called 'data'




##################################################
#   5. Fitting models using the Bayesian modeling software BUGS and JAGS
##################################################

### 5.1 Introduction

### 5.2 Introduction to BUGS software: WinBUGS, OpenBUGS, and JAGS

### 5.3 Linear model with normal response (normal GLM): multiple linear regression


#### Install first AHMbook package from CRAN
library(AHMbook)       # Load the package



# Generate data with data.fn from chapter 4
set.seed(24)
data <- data.fn()

#wait until data.fun() is finished
str(data)
attach(data)

# Summarize data by taking mean at each site and plot
(Cmean <- apply(data$C, 1, mean))
par(mfrow = c(1,3))
hist(Cmean, 50)               # Very skewed
plot(data$elev, Cmean)
plot(data$forest, Cmean)

# Package the data needed in a bundle
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = data$elev, forest = data$forest)
str(win.data)                    # Check what’s in win.data


# Write text file with model description in BUGS language
#when assigning priors, if fixing it at zero, same as commenting it out
cat(file = "multiple_linear_regression_model.txt",
"   
model {

# Priors - define for every parameter (either dnorm with large variance, or uniform dist); inverse in this case because precision is 1/variance
alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
tau <- pow(sd, -2)                   # Precision tau = 1/(sd^2)
sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale; have to define prior because define it in tau

# Likelihood
for (i in 1:M){
   Cmean[i] ~ dnorm(mu[i], tau)      # dispersion tau is precision (1/variance)
   mu[i] <- alpha0 + alpha1*elev[i] + alpha2*forest[i] + alpha3*elev[i]*forest[i]
}

# Derived quantities - residuals; could have defined under the likelihood for loop
for (i in 1:M){
   resi[i] <- Cmean[i] - mu[i]    
}
}"
)


# Initial values (have to give for at least some estimands)
inits <- function() list(alpha0 = rnorm(1,0,10), 
                         alpha1 = rnorm(1,0,10), 
                         alpha2 = rnorm(1,0,10), 
                         alpha3 = rnorm(1,0,10))


# Parameters monitored (i.e., for which estimates are saved)
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "resi")
#params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd")


# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3
# ni <- 10   ;   nt <- 1   ;   nb <- 0   ;  nc <- 8 # not run


# Call JAGS from R (ART <1 min)
library(jagsUI)
#?jags                 # Look at main function
out1 <- jags(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

# We can easily run JAGS on multiple cores by setting the argument parallel = TRUE:
out1 <- jags(win.data, inits, params, parallel = TRUE, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


par(mfrow = c(3,2))
traceplot(out1, param = c('alpha1', 'alpha2', 'resi[c(1,3, 5:6)]')) # Subset
traceplot(out1)                  # All params

# Overview of object created by jags()
names(out1)

# Summarize posteriors from JAGS run
print(out1, 2)
#mean is equivalent to MLE estimate
#sd is sd of MLE estimate?
#perc values are credible intervals
#f is proportion of MCMC samples that are on same side as the mean - one-sided test of significance
#Rhat - convergence
#DIC - bayes model selection; heuristic; shouldn't use if discrete random effects, which is most cases

str(out1)
#sims.list is actual output; 15,000 samples (3 chains * iterations): random sample from posterior distribution of intercept from linear regression
#sims.list$resi comes out as a matrix
#sims.list$deviance sorta like sums of squares

#all the rest of the outputs are summaries


# We can easily get various graphical overviews for both analyses (Fig. 5-2).
plot(out1)              # For JAGS analysis from jagsUI

#could grab top part of output and put in table in paper, or say "we report posterior mean distribution and sd"


par(mfrow = c(1, 2), mar = c(5,4,2,2), cex.main = 1)
whiskerplot(out1, param = c('alpha0', 'alpha1', 'alpha2', 'alpha3', 'sd')) #, 'resi[c(1,3, 5:7)]'))    # For JAGS analysis from jagsUI
library(denstrip)      # Similar, but more beautiful, with package denstrip
plot(out1$alpha0, xlim=c(-4, 4), ylim=c(1, 5), xlab="", ylab="", type="n", axes = F, main = "Density strip plots")
axis(1)
axis(2, at = 1:5, labels = c('alpha0','alpha1','alpha2','alpha3','sd'), las = 1)
abline(v = c(-4,-2,2,4), col = "grey")  ;  abline(v = 0)
for(k in 1:5){
   denstrip(unlist(out1$sims.list[k]), at = k, ticks = out1$summary[k, c(3,5,7)])
}

#similar to MLE estiamte for simple model such as this
(fm <- summary(lm(Cmean ~ elev*forest)))


print(cbind(out1$summary[1:5, 1:2], rbind(fm$coef[,1:2], c(fm$sigma, NA))), 4)
plot(lm(Cmean ~ elev*forest)) 
#looking at constant variance (homoscedastcity; would be rejected - shows higher variance when values are larger) and normal distr
#but shows should probably use Pois distr for counts rather than normal dist. on mean counts

#same plots using our JAGs output
mu <- out1$mean$alpha0 + out1$mean$alpha1 * elev + out1$mean$alpha2 * forest + out1$mean$alpha3 * elev * forest       # Compute the posterior mean of mu

par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
plot(1:M, out1$summary[6:272, 1], xlab = "Order of values", ylab = "Residual", frame.plot = F, ylim = c(-10, 15))
abline(h = 0, col = "red", lwd = 2)
segments(1:267, out1$summary[6:272, 3], 1:267, out1$summary[6:272, 7], col = "grey")
text(10, 14, "A", cex = 1.5)
hist(out1$summary[6:272, 1], xlab = "Residual", main = "", breaks = 50, col = "grey", xlim = c(-10, 15))
abline(v = 0, col = "red", lwd = 2)
text(-9, 48, "B", cex = 1.5)
qq <- qnorm(seq(0,0.9999,,data$M), mean = 0, sd = out1$summary[5, 1])
plot(sort(qq), sort(out1$summary[6:272, 1]), xlab = "Theoretical quantile", ylab = "Residual", frame.plot = F, ylim = c(-10, 15)) # could also use qqnorm()
abline(0, 1, col = "red", lwd = 2)
text(-4.5, 14, "C", cex = 1.5)
plot(mu, out1$summary[6:272, 1], xlab = "Predicted values", ylab = "Residual", frame.plot = F, ylim = c(-10, 15))
abline(h = 0, col = "red", lwd = 2)
segments(mu, out1$summary[6:272, 3], mu, out1$summary[6:272, 7], col = "grey")
text(-1, 14, "D", cex = 1.5)


confint(lm(Cmean ~ elev*forest))

#plotting the actual posterior samples for each covariate, in histograms with confidence interval
par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
hist(out1$sims.list$alpha1, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1$sims.list$alpha1, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(-2.4, 1.8, "A", cex = 1.5)
hist(out1$sims.list$alpha2, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1$sims.list$alpha2, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(1.7, 2, "B", cex = 1.5)
hist(out1$sims.list$alpha3, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1$sims.list$alpha3, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(-2.2, 1.2, "C", cex = 1.5)
hist(out1$sims.list$sd, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1$sims.list$sd, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(1.6, 4.9, "D", cex = 1.5)

library(coda)
HPDinterval(as.mcmc(out1$sims.list$sd), prob = 0.95)  # HPDI
quantile(out1$sims.list$sd, prob = c(0.025, 0.975))   # Percentile-based CRI

cbind(confint(lm(Cmean ~ elev*forest))[2:4,], out1$summary[2:4, c(3,7)])

#estimate of proportion that are greater than/less than 1.6
mean(out1$sims.list$alpha1 < -1.6)
mean(out1$sims.list$alpha1 < -1.6 & out1$sims.list$alpha1 > -1.8)

#not sure what this figure means - helps look at symmetry?
plot(out1$sims.list$alpha1, out1$sims.list$alpha2)
abline(h = c(2.5, 2.8), col = "red", lwd = 2)
abline(v = c(-1.9, -1.6), col = "red", lwd = 2)


mean(out1$sims.list$alpha1 < -1.6 & out1$sims.list$alpha1 > -1.9 & out1$sims.list$alpha2 > 2.5 & out1$sims.list$alpha2 < 2.8)

crazy.ratio <- out1$sims.list$alpha2 / abs(out1$sims.list$alpha1)
hist(crazy.ratio, main = "", breaks = 100, col = "grey", freq = F)
abline(v = quantile(crazy.ratio, prob = c(0.025, 0.975)), col = "red", lwd = 3)


mean(abs(out1$sims.list$alpha2 / out1$sims.list$alpha1) > 1)

# Compute expected abundance for a grid of elevation and forest cover
elev.pred <- seq(-1, 1,,100)                       # Values of elevation
forest.pred <- seq(-1,1,,100)                      # Values of forest cover
pred.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix
for(i in 1:100){
   for(j in 1:100){
      pred.matrix[i, j] <- out1$mean$alpha0 + out1$mean$alpha1 * elev.pred[i] + out1$mean$alpha2 * forest.pred[j] + out1$mean$alpha3 * elev.pred[i] * forest.pred[j]
   }
}

par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.main = 1.6, cex.axis = 1.5, cex.lab = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=elev.pred, y= forest.pred, z=pred.matrix, col = mapPalette(100), xlab = "Elevation", ylab = "Forest cover")
contour(x=elev.pred, y=forest.pred, z=pred.matrix, add = TRUE, lwd = 1, cex = 1.5)
title(main = "A")
matpoints(elev, forest, pch="+", cex=1.5)
abline(h = c(-1, -0.5, 0, 0.5, 1))

# Predictions for elev. at specific values of forest cover (-1,-0.5,0,0.5,1)
pred1 <- out1$mean$alpha0 + out1$mean$alpha1 * elev.pred + out1$mean$alpha2 * (-1) + out1$mean$alpha3 * elev.pred * (-1) 
pred2 <- out1$mean$alpha0 + out1$mean$alpha1 * elev.pred + out1$mean$alpha2 * (-0.5) + out1$mean$alpha3 * elev.pred * (-0.5)
pred3 <- out1$mean$alpha0 + out1$mean$alpha1 * elev.pred + out1$mean$alpha2 * 0 + out1$mean$alpha3 * elev.pred * 0
# pred3b <- out1$mean$alpha0 + out1$mean$alpha1 * elev.pred   # same
pred4 <- out1$mean$alpha0 + out1$mean$alpha1 * elev.pred + out1$mean$alpha2 * 0.5 + out1$mean$alpha3 * elev.pred * 0.5
pred5 <- out1$mean$alpha0 + out1$mean$alpha1 * elev.pred + out1$mean$alpha2 * 1 + out1$mean$alpha3 * elev.pred * 1
matplot(seq(-1, 1,,100), cbind(pred1, pred2, pred3, pred4, pred5), type = "l", lty= 1, col = "blue", ylab = "Prediction of mean count", xlab = "Elevation", ylim = c(-1.5, 7), lwd = 2)
title(main = "B")

#show degree of certainty for one of the predicted mean lines
#blue and red represent bayesian and frequentist predictions
pred.mat <- array(dim = c(length(elev.pred), length(out1$sims.list$alpha0)))
for(j in 1:length(out1$sims.list$alpha0)){
   pred.mat[,j] <- out1$sims.list$alpha0[j] + out1$sims.list$alpha1[j] * elev.pred + out1$sims.list$alpha2[j] * 0.5 + out1$sims.list$alpha3[j] * elev.pred * 0.5
}

CL <- apply(pred.mat, 1, function(x){quantile(x, prob = c(0.025, 0.975))})
plot(seq(-1, 1,,100), pred4, type = "l", lty= 1, col = "blue", ylab = "Prediction of mean count at forest = -0.5", xlab = "Elevation", las =1, ylim = c(-1.5, 7), lwd = 3)
matlines(seq(-1, 1,,100), t(CL), lty = 1, col = "blue", lwd = 2)
title(main = "C")

pred <- predict(lm(Cmean ~ elev*forest), newdata = data.frame(elev = seq(-1, 1,,100), forest = 0.5), se.fit = TRUE, interval = "confidence")
lines(seq(-1, 1,,100), pred$fit[,1], lty= 2, col = "red", lwd = 3)
matlines(seq(-1, 1,,100), pred$fit[,2:3], lty = 2, col = "red", lwd = 2)



### 5.6 Linear model with normal response (normal GLM): analysis of covariance (ANCOVA)
#(instead of forest being a continuous, bin into factor levels)

# Generate factor and plot raw data in boxplot as function of factor A
#in JAGs, factors must go from 1:N and be labeled as consecutive numbers
facFor <- as.numeric(forest < -0.5)         # Factor level 1
facFor[forest < 0 & forest > -0.5] <- 2     # Factor level 2
facFor[forest < 0.5 & forest > 0] <- 3      # Factor level 3
facFor[forest > 0.5] <- 4                   # Factor level 4
table(facFor)                               # every site assigned a level OK

#with the factors, Estimate Std.s are additive - gotta add together to get intercept for a given factor level
#e.g., intercept for factor 3 = 1+2+3
par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Cmean ~ factor(facFor), col = c("red", "blue", "green", "grey"), xlab = "Forest cover class", ylab = "Mean count of great tits", frame.plot = F, ylim = c(0,20))
text(0.8, 20, "A", cex=1.6)


# Bundle data
str(win.data <- list(Cmean = Cmean, M = length(Cmean), elev = elev, facFor = facFor))


# Specify model in BUGS language in effects parameterisation
cat(file = "ANCOVA1.txt","
model {

# Priors
alpha ~ dnorm(0, 1.0E-06)            # Prior for intercept = effect of level 1 of forest factor
beta2 ~ dnorm(0, 1.0E-06)            # Prior for slope = effect of elevation for level 1 of forest factor
beta1[1] <- 0                        # Set to zero effect of first level of facFor
beta3[1] <- 0                        # Set to zero effect of first level of facFor of elevation
for(k in 2:4){
   beta1[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
   beta3[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
}
tau <- pow(sd, -2)
sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

# Likelihood
for (i in 1:M){
   Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
  #alpha: global intercept
  #beta1: nested indexing: contribution of forest (vector of 4)
  #beta2: global slope
  #beta3: nested interaction factor
  #first four numbers describe intercepts, second four describe slopes; overparameterized as-is, so must fix intercept to 0
   mu[i] <- alpha + beta1[facFor[i]] + beta2 * elev[i] + beta3[facFor[i]] * elev[i] #four interaction factors
}
}
")

# Initial values
inits <- function() list(alpha = rnorm(1,,10), beta1 = c(NA, rnorm(3,,10)), beta2 = rnorm(1,,10), beta3 = c(NA, rnorm(3,,10)))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "sd")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R (ART <1 min)
out3 <- jags(win.data, inits, params, "ANCOVA1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(2, 2))
traceplot(out3)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev)))


# Summarize posteriors
print(out3, 3)

# Specify model in BUGS language
cat(file = "ANCOVA2.txt","
model {

# Priors
for(k in 1:4){
   alpha[k] ~ dnorm(0, 1.0E-06)       # Priors for intercepts
   beta[k] ~ dnorm(0, 1.0E-06)        # Priors for slopes
}
tau <- pow(sd, -2)
sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale

# Likelihood
for (i in 1:M){
   Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
   mu[i] <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
}

# Derived quantities: comparison of slopes (now you can forget the delta rule !)
for(k in 1:4){
   diff.vs1[k] <- beta[k] - beta[1]    # Differences relative to beta[1]
   diff.vs2[k] <- beta[k] - beta[2]    # ... relative to beta[2]
   diff.vs3[k] <- beta[k] - beta[3]    # ... relative to beta[3]
   diff.vs4[k] <- beta[k] - beta[4]    # ... relative to beta[4]
}
}
")

# Initial values
inits <- function() list(alpha = rnorm(4,,10), beta = rnorm(4,,10))

# Parameters monitored
params <- c("alpha", "beta", "sd", "diff.vs1", "diff.vs2", "diff.vs3", "diff.vs4")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R (ART <1 min) and summarize posteriors
out4 <- jags(win.data, inits, params, "ANCOVA2.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(2,2))   
traceplot(out4)
print(out4, 2)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev-1-elev)))

plot(elev[facFor==1], Cmean[facFor==1], col = "red", ylim = c(0, 20), xlab = "Elevation", ylab = "", frame.plot = F)
points(elev[facFor==2], Cmean[facFor==2], col = "blue")
points(elev[facFor==3], Cmean[facFor==3], col = "green")
points(elev[facFor==4], Cmean[facFor==4], col = "black")
abline(fm$coef[1,1], fm$coef[5,1], col = "red")
abline(fm$coef[2,1], fm$coef[6,1], col = "blue")
abline(fm$coef[3,1], fm$coef[7,1], col = "green")
abline(fm$coef[4,1], fm$coef[8,1], col = "black")
text(-0.8, 20, "B", cex=1.6)

# attach.bugs(out4)     # Allows to directly address the sims.list
str(out4$sims.list$diff.vs3)
par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(out4$sims.list$diff.vs3[,1], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.2, 0.8, "A", cex = 2)
hist(out4$sims.list$diff.vs3[,2], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.4, 0.8, "B", cex = 2)
hist(out4$sims.list$diff.vs3[,4], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-2.2, 0.8, "C", cex = 2)


# Prob. difference greater than 1
mean(out4$sims.list$diff.vs3[,1] > 1)
mean(out4$sims.list$diff.vs3[,2] > 1)
mean(out4$sims.list$diff.vs3[,4] > 1)
 [1] 0.6554667
 [1] 0.1981333
 [1] 0.003733333

 
 
 
 
### 5.8 Fitting a model with non-standard likelihood using the zeros or the ones tricks

# # Package the data needed in a bundle
# str(win.data <- list(Cmean1 = Cmean, Cmean2 = Cmean, zeros = rep(0, M), ones = rep(1, M), M = length(Cmean), elev = elev, forest = forest) ) # note 2 copies of response 
# 
# # Write text file with model description in BUGS language
# cat(file = "multiple_linear_regression_model.txt",
# "model {
# 
# # Priors
# for(k in 1:3){ # Loop over three ways to specify likelihood
#    alpha0[k] ~ dnorm(0, 1.0E-06)           # Prior for intercept
#    alpha1[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
#    alpha2[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
#    alpha3[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
#    sd[k] ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
# }
# var1 <- pow(sd[1], 2)                      # Variance in zeros trick
# var2 <- pow(sd[2], 2)                      # Variance in ones trick
# tau <- pow(sd[3], -2)                      # Precision tau = 1/(sd^2)
# 
# C1 <- 10000 # zeros trick: make large enough to ensure lam >= 0
# C2 <- 10000 # ones trick: make large enough to ensure p <= 1
# pi <- 3.1415926
# 
# # Three variants of specification of the likelihood
# for (i in 1:M){
# # 'Zeros trick' for normal likelihood
#    zeros[i] ~ dpois(phi[i])  # likelihood contribution is exp(-phi)
# #   negLL[i] <- log(sd[1]) + 0.5 * pow((Cmean1[i] - mu1[i]) / sd[1],2 )
#    negLL[i] <- -log(sqrt(1/(2*pi*var1))) + pow(Cmean1[i]-mu1[i],2)/(2*var1)
#    phi[i] <- negLL[i] + C1
#    mu1[i] <- alpha0[1] + alpha1[1]*elev[i] + alpha2[1]*forest[i] + alpha3[1]*elev[i]*forest[i]
# 
# # 'Ones trick' for normal likelihood
#    ones[i] ~ dbern(p[i])  # likelihood contribution is p directly
#    L[i] <- sqrt(1/(2*pi*var2)) * exp(-pow(Cmean1[i]-mu2[i],2)/(2*var2))
#    p[i] <- L[i] / C2
#    mu2[i] <- alpha0[2] + alpha1[2]*elev[i] + alpha2[2]*forest[i] + alpha3[2]*elev[i]*forest[i]
# 
# # Standard distribution function for the normal
#    Cmean2[i] ~ dnorm(mu3[i], tau)
#    mu3[i] <- alpha0[3] + alpha1[3]*elev[i] + alpha2[3]*forest[i] + alpha3[3]*elev[i]*forest[i]
# }
# }"
# )
# 
# # Initial values
# inits <- function() list(alpha0 = rnorm(3, 0, 10), alpha1 = rnorm(3,0,10), alpha2 = rnorm(3,0,10), alpha3 = rnorm(3,0,10))
# 
# # Parameters monitored (i.e., for which estimates are saved)
# params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd")
# 
# # MCMC settings
# ni <- 1200   ;   nt <- 1   ;   nb <- 200   ;  nc <- 3    # For JAGS
# 
# # Call JAGS
# library(jagsUI)
# outX <- jags(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
# print(outX)
# 
# 
# # Define negative log-likelihood function
# neglogLike <- function(param) {
#    alpha0 = param[1]
#    alpha1 = param[2]
#    alpha2 = param[3]
#    alpha3 = param[4]
#    sigma = exp(param[5])    # Estimate sigma on log-scale
#    mu = alpha0 + alpha1*elev + alpha2*forest + alpha3*elev*forest
# #   -sum(dnorm(Cmean, mean=mu, sd=sigma, log=TRUE))  # cheap quick way
# sum(-log(sqrt(1/(2*3.1415926*sigma^2))) + (Cmean-mu)^2/(2*sigma^2))
# }
# 
# # Find parameter values that minimize function value
# (fit <- optim(par = rep(0, 5), fn = neglogLike, method = "BFGS"))
#  [1]  1.6603052 -1.5764599  2.3440036 -0.8506793  0.6073662
# 
# exp(fit$par[5])            # Backtransform to get sigma




### 5.9 Poisson generalized linear model (Poisson GLM)

# Summarize data by taking max at each site
Cmax <- apply(C, 1, max)
table(Cmax)

# Bundle data
str(win.data <- list(Cmax = Cmax, M = length(Cmax), elev = elev, facFor = facFor, e = 0.0001) )

# Specify model in BUGS language
cat(file = "Poisson_GLM.txt","
model {

# Priors
for(k in 1:4){
   alpha[k] ~ dnorm(0, 1.0E-06)       # Prior for intercepts
   beta[k] ~ dnorm(0, 1.0E-06)        # Prior for slopes
}

# Likelihood; number of priors (above, 2*1:4 = 8) shows how many parameters
for (i in 1:M){
   Cmax[i] ~ dpois(lambda[i])         # note no variance parameter for Poisson
   log(lambda[i]) <- alpha[facFor[i]] +  beta[facFor[i]] * elev[i]  #log transformation = link function
   resi[i] <- (Cmax[i]-lambda[i]) / (sqrt(lambda[i])+e)   # Pearson residuals; standard error
} #number of parameters is 8 (no variance parameter w/ Pois)
}
")


# Initial values
inits <- function() list(alpha = rnorm(4,,3), beta = rnorm(4,,3))

# Parameters monitored
params <- c("alpha", "beta", "lambda", "resi")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R and summarize posteriors
out5 <- jags(win.data, inits, params, "Poisson_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(4,2))    ;    traceplot(out5, c("alpha[1:4]", "beta[1:4]"))
print(out5, 3)

par(mfrow = c(1, 3), mar = c(5,5,3,2), cex = 1.3, cex.lab = 1.5, cex.axis = 1.5)
hist(out5$summary[276:542, 1], xlab = "Pearson residuals", col = "grey", breaks = 50, main = "", freq = F, xlim = c(-5, 5), ylim = c(0, 0.57))
abline(v = 0, col = "red", lwd = 2)
text(-4.7, 0.54, "A", cex = 1.5)

plot(1:267, out5$summary[276:542, 1], main = "", xlab = "Order of data", ylab = "Pearson residual", frame.plot = F)
abline(h = 0, col = "red", lwd = 2)
text(8, 4, "B", cex = 1.5)

plot(out5$summary[9:275, 1],out5$summary[276:542, 1], main = "", xlab = "Predicted values", ylab = "Pearson residual", frame.plot = F, xlim = c(-1, 14))
abline(h = 0, col = "red", lwd = 2)
text(-0.5, 4, "C", cex = 1.5)

#same model in R
summary(glm(Cmax ~ factor(facFor)*elev-1-elev, family = poisson))

lambda2 <- array(dim = c(15000, 267))
for(j in 1:267){                            # Loop over sites
   lambda2[,j] <- exp(out5$sims.list$alpha[,facFor[j]] + out5$sims.list$beta[,facFor[j]] * elev[j]) # linear regression/backtransform
}
plot(out5$sims.list$lambda ~ lambda2, pch = ".")  # Check the two are identical
lm(c(out5$sims.list$lambda) ~ c(lambda2))

sorted.ele1 <- sort(elev[facFor == 1])
sorted.y1 <- out5$summary[9:275,][facFor == 1,][order(elev[facFor == 1]),]

# Plot A
par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(elev[facFor == 1], jitter(Cmax[facFor ==1]), ylab = "Maximum count", xlab = "Elevation (scaled)", frame.plot=F, ylim = c(0, 6))
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2) # Post. mean
lines(sorted.ele1, sorted.y1[,3], col = "grey", lwd = 2) # Lower 95% CL
lines(sorted.ele1, sorted.y1[,7], col = "grey", lwd = 2) # Upper 95% CL
text(-0.8, 6, "A", cex = 2)

# Plot B - shows variation with grey range
plot(sorted.ele1, sorted.y1[,1], type='n', xlab = "Elevation (scaled)", ylab = "", frame.plot = F, ylim = c(0, 6))
polygon(c(sorted.ele1, rev(sorted.ele1)), c(sorted.y1[,3], rev(sorted.y1[,7])), col='grey', border=NA)
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2)
text(-0.8, 6, "B", cex = 2)

# Plot C - shows variation - random sampling of slopes from posterior
elev.pred <- seq(-1,1, length.out = 200)  # Cov. for which to predict lambda
n.pred <- 50                             # Number of prediction profiles
pred.matrix <- array(NA, dim = c(length(elev.pred), n.pred))
for(j in 1:n.pred){
   sel <- sample(1:length(out5$sims.list$alpha[,1]),1) # Choose one post. draw
   pred.matrix[,j] <- exp(out5$sims.list$alpha[sel,1] + out5$sims.list$beta[sel,1] * elev.pred)
}
plot(sorted.ele1, sorted.y1[,1], type='n', xlab = "Elevation (scaled)", ylab = "", frame.plot = F, ylim = c(0, 6))
matlines(elev.pred, pred.matrix, col = "grey", lty = 1, lwd = 1)
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2)
text(-0.8, 6, "C", cex = 2)



### 5.11 Binomial generalised linear model (binomial GLM, logistic regression)

# Quantize counts from first survey and describe
y1 <- as.numeric(C[,1] > 0)  # Gets 1 if first count greater than zero
table(y1)

mean(N > 0)          # True occupancy
mean(y1)             # Observed occupancy after first survey


# Bundle data
str(win.data <- list(y1 = y1, M = length(y1), elev = elev, facFor = facFor) )

# Specify model in BUGS language
cat(file = "Bernoulli_GLM.txt","
model {

# Priors
for(k in 1:4){
   alpha[k] <- logit(mean.psi[k])     # intercepts
   mean.psi[k] ~ dunif(0,1)
   beta[k] ~ dnorm(0, 1.0E-06)        # slopes
}

# Likelihood
for (i in 1:M){
   y1[i] ~ dbern(theta[i])
   logit(theta[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
}
}
")

# Initial values
inits <- function() list(mean.psi = runif(4), beta = rnorm(4,,3))   # Priors 2

# Parameters monitored
params <- c("mean.psi", "alpha", "beta", "theta")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R (ART <1 min)
out6 <- jags(win.data, inits, params, "Bernoulli_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(4,2))    ;    traceplot(out6, c("alpha[1:4]", "beta[1:4]"))

print(out6, 2)

# Compare with MLEs
summary(glm(y1 ~ factor(facFor)*elev-1-elev, family = binomial))


# Plot of observed response vs. two covariates
par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
F1 <- facFor == 1 ; F2 <- facFor == 2 ; F3 <- facFor == 3 ; F4 <- facFor == 4
plot(jitter(y1,,0.05) ~ facFor, xlab = "Forest factor", ylab = "Observed occupancy probability", frame.plot = F, ylim = c(0, 1.15))
lines(1:4, out6$summary[1:4,1], lwd = 2)
segments(1:4, out6$summary[1:4,3], 1:4, out6$summary[1:4,7])
text(1.15, 1.1, "A", cex=1.6)

plot(elev[F1], jitter(y1,,0.1)[F1], xlab = "Elevation", ylab = "", col = "red", frame.plot = F)
points(elev[F2], jitter(y1,,0.05)[F2], col = "blue")
points(elev[F3], jitter(y1,,0.05)[F3], col = "green")
points(elev[F4], jitter(y1,,0.05)[F4], col = "grey")
lines(sort(elev[F1]), out6$mean$theta[F1][order(elev[F1])], col="red", lwd=2)
lines(sort(elev[F2]), out6$mean$theta[F2][order(elev[F2])], col="blue", lwd=2)
lines(sort(elev[F3]), out6$mean$theta[F3][order(elev[F3])], col="green", lwd=2)
lines(sort(elev[F4]), out6$mean$theta[F4][order(elev[F4])], col="grey", lwd=2)
text(-0.9, 1.1, "B", cex=1.6)



### 5.13 Random-effects Poisson GLM (Poisson GLMM)

# Bundle data
str(win.data <- list(C = data$C, M = nrow(data$C), J = ncol(data$C), elev = data$elev, 
                     forest = data$forest, elev.forest = data$elev * data$forest, wind = data$wind) )

# Specify model in BUGS language
cat(file = "RE.Poisson.txt","
model {

# Priors
mu.alpha ~ dnorm(0, 0.001)                # Mean hyperparam
tau.alpha <- pow(sd.alpha, -2)
sd.alpha ~ dunif(0, 10)                   # sd hyperparam
for(k in 1:4){
   alpha[k] ~ dunif(-10, 10)              # Regression params
}

# Likelihood
#intercept, slopes for elev and wind speed; but counts aren't independent of site (random site effect)
#random site effect - directly index intercept (alpha0[i]) = fixed, but assign distribution to make it random
#thus the line in the prior line alpha0[i] ~ dnorm(mu, tau)
#estimates different for each combination of i and j [i,j]; lambda is mean central structure
#alpha0 is the only one indexed to be different for each site and occassion
for (i in 1:M){
   alpha0[i] ~ dnorm(mu.alpha, tau.alpha) # Random effects and hyperparams
   re0[i] <- alpha0[i] - mu.alpha         # zero-centered random effects
   for(j in 1:J){
      C[i,j] ~ dpois(lambda[i,j]) 
      log(lambda[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] + 
alpha[3] * elev.forest[i] + alpha[4] * wind[i,j]
   }
}
}")

inits <- function() list(alpha0 = rnorm(M), alpha = rnorm(4)) # Inits
params <- c("mu.alpha", "sd.alpha", "alpha0", "alpha", "re0") # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call JAGS from R (ART 6-7 min) and summarize posteriors
out8 <- jags(win.data, 
             #inits, 
             parameters.to.save = params, model.file = "RE.Poisson.txt", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow = c(3,2))  ;  traceplot(out8, c("mu.alpha", "sd.alpha", "alpha[1:3]"))

print(out8, 3)

Cvec <- as.vector(data$C)            # Vector of M*J counts
elev.vec <- rep(data$elev, data$J)        # Vectorized elevation covariate
forest.vec <- rep(data$forest, data$J)    # Vectorized forest covariate
wind.vec <- as.vector(data$wind)     # Vectorized wind covariate
fac.site <- factor(rep(1:data$M, data$J)) # Site indicator (factor)
cbind(Cvec, fac.site, elev.vec, forest.vec, wind.vec) # Look at data

# Fit same model using maximum likelihood (NOTE: glmer uses ML instead of REML)
library(lme4)
summary(fm <- glmer(Cvec ~ elev.vec*forest.vec + wind.vec + (1| fac.site), family = poisson))              # Fit model
ranef(fm)                       # Print zero-centered random effects

# We compare the fixed-effects estimates in a table below, and also the random-effects estimates from the Bayesian and the non-Bayesian analyses in a table and a graph (Fig. 5–17A).

# Compare fixed-effects estimates, Bayesian post. and freq LMEs
print(cbind(out8$summary[c(1:2, 270:273), 1:2], rbind(summary(fm)$coef[1,1:2], c(sqrt(summary(fm)$varcor$fac.site), NA), summary(fm)$coef[c(2,3,5,4),1:2])), 3)

# Compare graphically non-Bayesian and Bayesian random effects estimates
Freq.re <- ranef(fm)$fac.site[,1]         # Non-Bayesian estimates (MLEs)
Bayes.re <- out8$summary[274:540,]        # Bayesian estimates

par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Freq.re, Bayes.re[,1], xlab = "Non-Bayesian (glmer)", ylab = "Bayesian (BUGS)", xlim = c(-0.4, 0.42), ylim = c(-2, 2), frame.plot = F, type = "n")
segments(Freq.re, Bayes.re[,3], Freq.re, Bayes.re[,7], col = "grey", lwd = 0.5)
abline(0, 1, lwd = 2)
points(Freq.re, Bayes.re[,1])
text(-0.38, 2, "A", cex=1.6)


wind.pred <- seq(-1, 1, , 1000)     # Covariate values for prediction
pred <- array(NA, dim = c(1000, 267))
for(i in 1:267){
   pred[,i]<- exp(out8$mean$alpha0[i] + out8$mean$alpha[1] * 0 + out8$mean$alpha[2] * 0 + out8$mean$alpha[3] * 0 + out8$mean$alpha[4] * wind.pred)    # Predictions for each site
}

matplot(wind.pred, pred, type = "l", lty = 1, col = "grey", xlab = "Wind speed", ylab = "Expected count", frame.plot = F, ylim = c(0, 4))
lines(wind.pred, exp(out8$mean$mu.alpha + out8$mean$alpha[4] * wind.pred), col = "black", lwd = 3)
text(-0.9, 4, "B", cex=1.6)



### 5.14 Random-effects binomial GLM (binomial GLMM)

# Get detection/nondetection response
y <- C
y[y > 0] <- 1


# Bundle data
win.data <- list(y = y, M = nrow(y), J = ncol(y), elev = elev, forest = forest, elev.forest = elev * forest, wind = wind)
str(win.data)

# Specify model in BUGS language
cat(file = "RE.Bernoulli.txt","
model {

# Priors
mu.alpha0 <- logit(mean.theta)              # Random intercepts
mean.theta ~ dunif(0,1)
tau.alpha0 <- pow(sd.alpha0, -2)
sd.alpha0 ~ dunif(0, 10)
mu.alpha4 ~ dnorm(0, 0.001)                 # Random slope on wind
tau.alpha4 <- pow(sd.alpha4, -2)
sd.alpha4 ~ dunif(0, 10)
for(k in 1:3){
   alpha[k] ~ dnorm(0, 0.001)               # Slopes
}

# Likelihood
for (i in 1:M){
   alpha0[i] ~ dnorm(mu.alpha0, tau.alpha0) # Intercept random effects
   re00[i] <- alpha0[i] - mu.alpha0         # same zero-centered
   alpha4[i] ~ dnorm(mu.alpha4, tau.alpha4) # Slope random effects
   re04[i] <- alpha4[i] - mu.alpha4         # same zero-centered
   for(j in 1:J){
      y[i,j] ~ dbern(theta[i,j])
      logit(theta[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] + alpha[3] * elev.forest[i] + alpha4[i] * wind[i,j]
   }
}
}")

# Other model run preparations
inits <- function() list(alpha0 = rnorm(M), alpha4 = rnorm(M))# Inits
params <- c("mu.alpha0", "sd.alpha0", "alpha0", "alpha", "mu.alpha4", "sd.alpha4", "alpha4", "re00", "re04")                        # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call JAGS from R (ART 2.5 min)
out9 <- jags(win.data, inits, params, "RE.Bernoulli.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(2,2))
traceplot(out9, c("mu.alpha0", "sd.alpha0", "alpha[1:3]", "mu.alpha4", "sd.alpha4"))
print(out9, 3)

yvec <- as.vector(y)            # Vector of M*J counts
elev.vec <- rep(elev, J)        # Vectorized elevation covariate
forest.vec <- rep(forest, J)    # Vectorized forest covariate
wind.vec <- as.vector(wind)     # Vectorized wind covariate
fac.site <- factor(rep(1:M, J)) # Site indicator (factor)
cbind(yvec, fac.site, elev.vec, forest.vec, wind.vec) # Look at data

# Fit same model using maximum likelihood
library(lme4)                   # Load package
summary(frem <- glmer(yvec ~ elev.vec*forest.vec + wind.vec + (wind.vec || fac.site), family = binomial))              # Fit model


# Compare Bayesian and non-Bayesian estimates
print(out9$summary[c(1:2, 270:274),c(1:3,7:9)], 4)

(re <- ranef(frem))                 # Print zero-centered random effects