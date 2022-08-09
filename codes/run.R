## Load function for simulating data
## sim3A.R, sim3B.R, or sim3C.R
source("sim3C.R")

## Load function for ensemble 
source("ranger-addon.R")

## Load packages
library(ranger)
library(survival)

## Generate training data from scenario 3C
dat <- simDat3C(400, .2)

## Fits ranger
fit <- ranger(Surv(Time, status) ~ ., data = dat)

## Generate testing data and computes the predicted probabilities
dat0 <- simDat3C(500, 0)
predS <- getSurv(fit, dat, dat0)
trueS <- getTrueSurv3C(dat0)

t0 <- seq(0, quantile(dat0$Time, .9), .01)
## Intergrated absolute error
mean(sapply(1:nrow(dat0), function(i) abs(predS[[i]](t0) - trueS[[i]](t0))))
## Intergrated MSE
mean(sapply(1:nrow(dat0), function(i) (predS[[i]](t0) - trueS[[i]](t0))^2))


## Variable importance

vnames <- setdiff(names(dat), c("Time", "status", "D"))
p <- length(vnames)

getVi <- function(fit, df) {
  df <- df[order(df$Time + df$status),]
  tt <- sort(df$Time)
  rownames(df) <- NULL
  sapply(1:nrow(df), function(s) {
    fit.tmp <- fit
    takei <- which(unlist(lapply(fit$forest$OobSampleIDs, function(x) any(x %in% (s - 1)))))
    fit.tmp$forest$OobSampleIDs <- fit.tmp$forest$OobSampleIDs[takei]
    fit.tmp$forest$num.trees <- fit.tmp$num.trees <- length(takei)
    fit.tmp$forest$child.nodeIDs <- fit.tmp$forest$child.nodeIDs[takei]
    fit.tmp$forest$split.values <- fit.tmp$forest$split.values[takei]
    fit.tmp$forest$SampleIDs <- fit.tmp$forest$SampleIDs[takei]
    fit.tmp$forest$chf <- fit.tmp$forest$chf[takei]
    fit.tmp$forest$split.varIDs <- fit.tmp$forest$split.varIDs[takei]
    getSurv(fit.tmp, df, df[s,])
  })
}


sc <- with(survfit(Surv(Time, 1 - status) ~ 1, data = dat), stepfun(time, c(1, surv)))
getCON <- function(tt, si, sc, df0) {
  c1 <- outer(df0$status, rep(1, length(df0$status)))
  c2 <- outer(df0$Time <= tt, df0$Time > tt)
  l <- length(si)
  stmp <- sapply(1:l, function(x) si[[x]](tt))
  c3 <- 1 * outer(stmp, stmp, "<") + 0.5 * outer(stmp, stmp, "==")
  c4 <- outer(sc(df0$Time + df0$D), sc(df0$Time + tt))
  c4 <- ifelse(c4 < 1e-5, NA, c4)
  c4 <- 1
  sum(c1 * c2 * c3 / c4, na.rm = TRUE) / sum(c1 * c2 / c4, na.rm = TRUE)
}

## Based on testing data
mean(sapply(t0, getCON, predS, sc, dat0), na.rm = T)

## Based on training data
(c0 <- mean(sapply(t0, getCON, getSurv(fit, dat, dat), sc, dat), na.rm = T))

i <- 2
B <- 5

now <- Sys.time()
replicate(B, {
  dat2 <- dat[order(dat$Time),]
  dat2[,vnames[i]] <- sample(dat2[,vnames[i]])
  f <- getVi(fit, dat2)
  mean(sapply(t0, getCON, f, dat2), na.rm = T)})
Sys.time() - now


vimps <- sapply(1:length(vnames), function(i) 
  replicate(B, {
    dat2 <- dat[order(dat$Time),]
    dat2[,vnames[i]] <- ifelse(dat$status > 0, sample(dat2[,vnames[i]]), dat2[,vnames[i]])
    f <- getVi(fit, dat2)
    mean(sapply(t0, getCON, f, sc, dat2), na.rm = T)}))

library(xtable)
library(dplyr)
library(forcats)
library(ggridges)
library(ggpubr)

dd <- data.frame(vars = rep(vnames, each = B), vimp = c0 - c(vimps))
dd <- dd %>% mutate(vars = fct_reorder(vars, vimp, .fun = 'median'))

dd2 <- do.call(rbind, lapply(split(dd, dd$vars), function(d)
  data.frame(vars = unique(d$vars), vimp = median(d$vimp))))
ggplot(dd2, aes(x = vars, y = vimp)) + geom_bar(stat = "identity") + coord_flip()

## ridge
ggplot(dd, aes(x = vimp, y = vars)) + geom_density_ridges(rel_min_height = 1e-3)
  ## geom_density_ridges(bandwidth = 1e-2, alpha = .5, rel_min_height = 5e-3) 

## boxplot
ggplot(dd, aes(x = vars, y = vimp)) + geom_boxplot() + coord_flip()

library(mcreplicate)

vimps <- sapply(1:length(vnames), function(i) 
  mc_replicate(B, {
    dat2 <- dat[order(dat$Time),]
    dat2[,vnames[i]] <- ifelse(dat$status > 0, sample(dat2[,vnames[i]]), dat2[,vnames[i]])
    f <- getVi(fit, dat2)
    mean(sapply(t0, getCON, f, sc, dat2), na.rm = T)}))

dd <- data.frame(vars = rep(vnames, each = B), vimp = con0 - c(vimps))
dd %>% mutate(vars = fct_reorder(vars, vimp, .fun = 'median')) %>%
  ggplot(aes(x = vars, y = vimp)) + geom_boxplot() + coord_flip()




vimps <- sapply(1:2, function(i) 
  mc_replicate(B, {
    dat2 <- dat[order(dat$Time),]
    dat2[,vnames[i]] <- ifelse(dat$status > 0, sample(dat2[,vnames[i]]), dat2[,vnames[i]])
    f <- getVi(fit, dat2)
    mean(sapply(t0, getCON, f, sc, dat2), na.rm = T)}))
