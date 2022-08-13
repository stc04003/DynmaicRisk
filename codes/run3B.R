## Load function for simulating data
source("sim3B.R")

## Load function for ensemble 
source("ranger-addon.R")

## Load packages
library(ranger)
library(survival)

## Generate training data from scenario 3C
set.seed(0); dat <- simDat3B(n = 400, cen = .2)
tmp <- dat[,colSums(is.na(dat)) > 0]
tmp[is.na(tmp)] <- 1e5
dat[is.na(dat)] <- -1e5
dat <- data.frame(dat, tmp)
names(dat)[grep("W", names(dat))] <- 
  paste0(paste0("W.", 1:10), ".", rep(1:2, each = 10))

## Fits ranger
fit <- ranger(Surv(Time, status) ~ ., data = dat, min.node.size = 15)

## Generate testing data and computes the predicted probabilities
dat0 <- simDat3B(500, 0)
tmp <- dat0[, colSums(is.na(dat0)) > 0]
tmp[is.na(tmp)] <- 1e5
dat0[is.na(dat0)] <- -1e5
dat0 <- data.frame(dat0, tmp)
names(dat0)[grep("W", names(dat0))] <- 
  paste0(paste0("W.", 1:10), ".", rep(1:2, each = 10))

predS <- getSurv(fit, dat, dat0)
head(predS, 3)
trueS <- getTrueSurv3B(dat0)

t0 <- seq(0, quantile(dat0$Time[dat0$status > 0], .9), .01)
## Intergrated absolute error
mean(sapply(1:nrow(dat0), function(i) abs(predS[[i]](t0) - trueS[[i]](t0))))
## Intergrated MSE
mean(sapply(1:nrow(dat0), function(i) (predS[[i]](t0) - trueS[[i]](t0))^2))

## Intergrated CON
sc <- with(survfit(Surv(Time + D, 1 - status) ~ 1, data = dat, subset = Time < 1),
           stepfun(time, c(1, surv)))
mean(sapply(t0, getCON, predS, sc, dat0), na.rm = TRUE)

## Intergrated CON for the training data
f <- getSurv(fit, dat, dat)
(con0 <- mean(sapply(t0, getCON, f, sc, dat), na.rm = T))
f <- getSi(fit, dat)
(con0 <- mean(sapply(t0, getCON, f, sc, dat), na.rm = T))

vnames <- c(paste("Z", 1:10, sep = "."), paste0("W.", 1:10, "."))
vnames <- paste("Z", 1:10, sep = ".")
B <- 5

## Variable importance 
set.seed(0)
library(mcreplicate)
vimps <- sapply(1:length(vnames), function(i) 
  mc_replicate(B, {
    dat2 <- dat
    if (substr(vnames[i], 1, 1) == "Z")
      dat2[,vnames[i]] <- sample(dat2[,vnames[i]])
    else {
      permID <- grep(vnames[i], names(dat2))
      toPerm <- dat2[,permID[1]] > -1e5
      dat2[toPerm, permID] <- dat2[sample(which(toPerm)), permID]
    }
    f <- getSi(fit, dat2)
    mean(sapply(t0, getCON, f, sc, dat2), na.rm = T)
}))

## Variable importance plot
library(dplyr)
library(forcats)
library(ggplot2)

dd <- data.frame(vars = rep(vnames, each = B), vimp = con0 - c(vimps))
dd %>% mutate(vars = fct_reorder(vars, vimp, .fun = 'median')) %>% 
  ggplot(aes(x = vars, y = vimp)) + geom_boxplot() + coord_flip()


dd <- data.frame(vars = vnames, vimp = con0 - colMeans(vimps))
dd %>% mutate(vars = fct_reorder(vars, vimp)) %>% 
  ggplot(aes(x = vars, y = vimp)) + geom_bar(stat = "identity") + coord_flip()


vnames <- paste("Z", 1:10, sep = ".")
B <- 15

f <- getSi(fit, dat)
(con0 <- mean(sapply(t0, getCON2, f, sc, dat), na.rm = T))

## Variable importance for Z
set.seed(0)
vimps0 <- sapply(1:length(vnames), function(i) 
  mc_replicate(B, {
    dat2 <- dat ## [order(dat$Time),]
    permID <- grep(vnames[i], names(dat2))[1]
    dat2[,permID] <- sample(dat2[,permID])
    ## toPerm <- dat2[,permID[1]] > -1e5
    ## dat2[toPerm, permID] <- dat2[sample(which(toPerm)),permID]
    f <- getSi(fit, dat2)
    mean(sapply(t0, getCON2, f, sc, dat2), na.rm = T)
  }))

## boxplot
dd <- data.frame(vars = rep(vnames, each = B), vimp = con0 - c(vimps0))
dd %>% mutate(vars = fct_reorder(vars, vimp, .fun = 'median')) %>% 
  ggplot(aes(x = vars, y = vimp)) + geom_boxplot() + coord_flip()

## barplot
con0 <- mean(sapply(t0, getCON2, getSi(fit, dat), sc, dat), na.rm = T)
dd <- data.frame(vars = vnames, vimp = con0 - colMeans(vimps0))
dd %>% mutate(vars = fct_reorder(vars, vimp, .fun = 'median')) %>% 
  ggplot(aes(x = vars, y = vimp)) + geom_bar(stat = "identity") + coord_flip()



f <- getSi(fit, dat)

sc <- with(survfit(Surv(runif(1e5, 0, 5), rep(1, 1e5)) ~ 1), stepfun(time, c(1, surv)))
(con0 <- mean(sapply(t0, getCON, f, sc, dat), na.rm = T))
