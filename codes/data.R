#; Function used to generate W(t)
Wt <- function(t, a) a * pweibull(.2 * t, 2)

#' Function to generate simulated data
#'
#' Generate data from irreversible multi-state models with three states: healthy, diseased, and death.
#' The landmark time is assumed at the intermediate event.
#'
#' 
simDat <- function(n) {
  p <- 10
  Z <- matrix(rnorm(n * p), n)
  a <- matrix(runif(n * p, -1, 1), n)
  e <- rnorm(n)
  g <- rgamma(n, 2 , 2)
  t1 <- runif(n, 0, 5)
  t2 <- exp(-5 + rowSums(Wt(t1, a[,1:3])) + rowSums(Z[,1:3] * Z[,1:3]) + rowSums(Wt(t1, a[,1:3]) * Z[,1:3]) + log(1 + t1) + g + e)
  pb <- 1 / (1 + 1 / exp(rowSums(Wt(t1, a[,1:3])) + rowSums(Z[,1:3]) + g))
  dat0 <- data.frame(id = 1:n, Time = t1 + t2, baseCov, W = Wt(t1, a)) 
  dat <- dat0[rbinom(n, 1, pb) > 0,]
  dat$id <- 1:nrow(dat)
  dat
}
