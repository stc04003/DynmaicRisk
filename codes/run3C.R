## #####################################################################
## R codes used in ../Examples/sim3C.Rmd
## #####################################################################

source("../codes/sim3C.R")
set.seed(1); dat <- simDat3C(n = 400, cen = .2)

library(survival)
library(ranger)
(fit <- ranger(Surv(Time, status) ~ ., data = dat,
               keep.inbag = TRUE, min.node.size = 15))

source("../codes/ranger-addon.R")
dat0 <- simDat3C(500, 0)
predS <- getSurv(fit, dat, dat0)
head(predS, 3)
trueS <- getTrueSurv3C(dat0)

t0 <- seq(0, quantile(dat0$Time[dat0$status > 0], .9), .01)
## Intergrated absolute error
mean(sapply(1:nrow(dat0), function(i) abs(predS[[i]](t0) - trueS[[i]](t0))))
## Intergrated MSE
mean(sapply(1:nrow(dat0), function(i) (predS[[i]](t0) - trueS[[i]](t0))^2))
sc <- with(survfit(Surv(Time + D, 1 - status) ~ 1, data = dat), 
           stepfun(time, c(1, surv)))
## Intergrated CON
mean(sapply(t0, getCON, predS, sc, dat0), na.rm = TRUE)

con0 <- mean(sapply(t0, getCON2, getSi(fit, dat), sc, dat), na.rm = T)
vnames <- setdiff(names(dat), c("Time", "status", "D"))
B <- 100
library(mcreplicate)
set.seed(1)
vimps <- sapply(1:length(vnames), function(i) 
  mc_replicate(B, {
    dat2 <- dat
    dat2[,vnames[i]] <- sample(dat2[,vnames[i]])
    f <- getSi(fit, dat2)
    mean(sapply(t0, getCON2, f, sc, dat2), na.rm = T)}))

library(dplyr)
library(forcats)
library(ggplot2)
dd <- data.frame(vars = rep(vnames, each = B), vimp = con0 - c(vimps))
dd %>% mutate(vars = fct_reorder(vars, vimp, .fun = 'median')) %>% 
  ggplot(aes(x = vars, y = vimp)) + geom_boxplot() + 
  coord_flip() + xlab("") + ylab("")
