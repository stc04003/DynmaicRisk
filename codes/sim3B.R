#' Function used to generate W(t)
Wt <- function(t, a) a * pweibull(.2 * t, 2)

#' Function to generate simulated data
#'
#' Generate data from irreversible multi-state models with three states: healthy, diseased, and death.
#' The landmark time is 2 (administrative censoring at 3)
#' and longitudinal variables observed at the intermediate event
#'
#' @param n number of subjects at baseline
#' @param cen is the censoring percentage at baseline; available values are 0%, 20%, and 40%.
simDat3B <- function(n, cen) {
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
  dat0 <- data.frame(Time = t1 + t2 * (rbinom(n, 1, pb) > 0),
                     status = 1, D = t1, Z = Z, W = Wt(t1, a))
  dat0$Time <- pmin(dat0$Time, cc, 3)
  dat0$status <- 1 * (dat0$Time < pmin(cc, 3))
  dat0$D <- pmin(dat0$D, 2)
  dat <- subset(dat0, Time > 2)
  dat$Time <- dat$Time - 2
  dat[dat$D == 2, grep("W", names(dat0))] <- NA
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
getTrueSurv3B <- function(dat) {
  trSurv <- list()
  for (i in 1:nrow(dat)) {
    n2 <- 1e5
    W0 <- as.numeric(dat[i, grep("W", names(dat))])[1:3]
    Z0 <- as.numeric(dat[i, grep("Z", names(dat))])[1:3]    
    g <- rgamma(n2, 2, 2)
    pb2 <- 1 / (1 + 1 / exp(sum(W0) + sum(Z0) + g))
    pi2 <- rbinom(n2, 1, pb2)
    t1 <- as.numeric(dat$D[i])
    if (t1 == 2) {
      t12 <- runif(n2, 0, 5)
      t22 <- ifelse(pi2 == 0, t12, 
                    t12 + exp(-5 + sum(W0) + sum(Z0^2) + sum(W0 * Z0) + log(1 + t12) + g + rnorm(n2)))
      t22 <- sort(t22[t12 > 2] - 2)
    } else {
      t22 <- t1 + exp(-5 + sum(W0) + sum(Z0^2) + sum(W0 * Z0) + log(1 + t1) + g + rnorm(n2))
      t22 <- sort(t22[pi2 > 0] - 2)      
    }
    trSurv[[i]] <- stepfun(t22, c(1, 1 - ecdf(t22)(t22)))      
  }
  trSurv
}
