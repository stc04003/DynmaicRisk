#' Function used to generate W(t)
Wt <- function(t, a) a * pweibull(.2 * t, 2)

#' Function to generate simulated data
#'
#' Generate data from irreversible multi-state models with three states: healthy, diseased, and death.
#' The landmark time is 2 (administrative censoring at 3)
#' and longitudinal variables observed at times 1, 2, 3, 4, and 5
#'
#' @param n number of subjects at baseline.
#' @param cen is the censoring percentage at baseline; available values are 0%, 20%, and 40%.
#' @param test indicates whether we want to include "a" in output.
simDat3A <- function(n, cen, test = FALSE) {
  p <- 10
  Z <- matrix(rnorm(n * p), n)
  a <- matrix(runif(n * p, -1, 1), n)
  e <- rnorm(n)
  g <- rgamma(n, 2 , 2)
  t1 <- runif(n, 0, 5)
  t2 <- exp(-5 + rowSums(Wt(t1, a[,1:3])) + rowSums(Z[,1:3] * Z[,1:3]) +
            rowSums(Wt(t1, a[,1:3]) * Z[,1:3]) + log(1 + t1) + g + e)
  pb <- 1 / (1 + 1 / exp(rowSums(Wt(t1, a[,1:3])) + rowSums(Z[,1:3]) + g))
  if (cen == 0) cc <- rep(Inf, n)
  if (cen == .2) cc <- rexp(n, .04)
  if (cen == .4) cc <- rexp(n, .12)
  if (test)
    dat0 <- data.frame(Time = t1 + t2, status = 1, D = t1, Z = Z, a = a, 
                       W1 = Wt(1, a), W2 = Wt(2, a), W3 = Wt(3, a), W4 = Wt(4, a), W5 = Wt(5, a))
  else
    dat0 <- data.frame(Time = t1 + t2, status = 1, D = t1, Z = Z,
                       W1 = Wt(1, a), W2 = Wt(2, a), W3 = Wt(3, a), W4 = Wt(4, a), W5 = Wt(5, a))
  keep <- rbinom(n, 1, pb) > 0
  dat0$Time <- pmin(dat0$Time, cc)
  dat0$status <- 1 * (dat0$Time < cc)
  dat <- subset(dat0[keep,], Time >= D)
  for (i in 1:nrow(dat)) {
    if (floor(dat$Time[i]) < 5) dat[i, grep("W5", names(dat0))] <- NA
    if (floor(dat$Time[i]) < 4) dat[i, grep("W4", names(dat0))] <- NA
    if (floor(dat$Time[i]) < 3) dat[i, grep("W3", names(dat0))] <- NA
    if (floor(dat$Time[i]) < 2) dat[i, grep("W2", names(dat0))] <- NA
    if (floor(dat$Time[i]) < 1) dat[i, grep("W1", names(dat0))] <- NA
  }
  dat$Time <- dat$Time - dat$D
  dat <- dat[order(dat$Time),]
  rownames(dat) <- NULL
  return(dat)
}


#' Function to calculate (empirically) true landmark probability
#' given the dataset is generated from the simDat() function above
#'
#' @param dat dataset generated from simDat()
#'
#' This returns a list; the ith element in the list is the
#' landmark probability for the ith row of data
getTrueSurv3A <- function(dat) {
  trSurv <- list()
  for (i in 1:nrow(dat)) {
    n2 <- 1e5
    a0 <- as.numeric(dat[i, grep("a", names(dat))])[1:3]
    W0 <- Wt(dat$D[i], a0)
    Z0 <- as.numeric(dat[i, grep("Z", names(dat))])[1:3]    
    t1 <- as.numeric(dat$D[i])
    g <- rgamma(n2, 2, 2)
    pb2 <- 1 / (1 + 1 / exp(sum(W0) + sum(Z0) + g))
    pi2 <- rbinom(n2, 1, pb2)
    t22 <- exp(-5 + sum(W0) + sum(Z0^2) + sum(W0 * Z0) + log(1 + t1) + g + rnorm(n2))
    t22 <- sort(t22[pi2 > 0])
    trSurv[[i]] <- stepfun(t22, c(1, 1 - ecdf(t22)(t22)))
  }
  trSurv
}
