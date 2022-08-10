## Load function for simulating data
source("sim3C.R")

## Load function for ensemble 
source("ranger-addon.R")

## Load packages
library(ranger)
library(survival)

## Generate training data from scenario 3C
set.seed(1); dat <- simDat3C(n = 400, cen = .2)

## Fits ranger
fit <- ranger(Surv(Time, status) ~ ., data = dat)

## Generate testing data and computes the predicted probabilities
dat0 <- simDat3C(500, 0)
predS <- getSurv(fit, dat, dat0)
head(predS, 3)
trueS <- getTrueSurv3C(dat0)

t0 <- seq(0, quantile(dat0$Time[dat0$status > 0], .9), .01)
## Intergrated absolute error
mean(sapply(1:nrow(dat0), function(i) abs(predS[[i]](t0) - trueS[[i]](t0))))
## Intergrated MSE
mean(sapply(1:nrow(dat0), function(i) (predS[[i]](t0) - trueS[[i]](t0))^2))

sc <- with(survfit(Surv(Time + D, 1 - status) ~ 1, data = dat), stepfun(time, c(1, surv)))
mean(sapply(t0, getCON, predS, sc, dat0), na.rm = TRUE)

## Intergrated CON for the training data
(con0 <- mean(sapply(t0, getCON, getSurv(fit, dat, dat), sc, dat), na.rm = T))

vnames <- setdiff(names(dat), c("Time", "status", "D"))
B <- 100

## Variable importance 
set.seed(0)
library(mcreplicate)
vimps <- sapply(1:length(vnames), function(i) 
  mc_replicate(B, {
    dat2 <- dat[order(dat$Time),]
    dat2[,vnames[i]] <- ifelse(dat2$status > 0, sample(dat2[,vnames[i]]), dat2[,vnames[i]])
    f <- getSi(fit, dat2)
    mean(sapply(t0, getCON, f, sc, dat2), na.rm = T)}))

## Variable importance plot
library(dplyr)
library(forcats)
library(ggplot2)
dd <- data.frame(vars = rep(vnames, each = B), vimp = con0 - c(vimps))
dd %>% mutate(vars = fct_reorder(vars, vimp, .fun = 'median')) %>% 
  ggplot(aes(x = vars, y = vimp)) + geom_boxplot() + coord_flip()
