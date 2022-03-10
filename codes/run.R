## Load function for simulating data
source("data.R")

## Load function for ensemble 
source("ranger-addon.R")

## Load packages
library(ranger)
library(survival)

## Generate training data and fitting
dat <- simDat(400, .2)
fit <- ranger(Surv(Time, status) ~ ., data = dat)

## Generate testing data and computes the predicted probabilities
dat0 <- simDat(500, 0)
predS <- getSurv(fit, dat, dat0)
trueS <- getTrueSurv(dat0)
  
t0 <- seq(0, quantile(dat0$Time, .9), .01)
## Intergrated absolute error
mean(sapply(1:nrow(dat0), function(i) abs(predS[[i]](t0) - trueS[[i]](t0))))
## Intergrated MSE
mean(sapply(1:nrow(dat0), function(i) (predS[[i]](t0) - trueS[[i]](t0))^2))
