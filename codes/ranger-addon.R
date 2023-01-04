#' Function used to compute the predicted survival using ensemble method
#'
#' @param fit is a ranger fit 
#' @param trainDat is the training data used to create ranger
#' @param testDat is the testing data to computed the predicted survival
#'
#' This returns a list; the ith element in the list is the
#' landmark probability for the ith row of the testing data
getSurv <- function(fit, trainDat, testDat) {
  resChar <- strsplit(as.character(fit$call)[2], "\\s~")[[1]][1] 
  resChar <- gsub("Surv|[(]|[)]", "", resChar)
  resChar <- strsplit(resChar, ", ")[[1]]
  Time <- as.numeric(unlist(trainDat[resChar[1]]))
  Status <- as.numeric(unlist(trainDat[resChar[2]]))
  n <- nrow(trainDat)
  wbij <- do.call(cbind, fit$inbag.counts)
  trainDatNodes <- predict(fit, trainDat, type = "terminalNodes")$predictions
  testDatNodes <- predict(fit, testDat, type = "terminalNodes")$predictions
  predSV <- sapply(1:nrow(testDat), function(i) {
    wij <- rowSums(wbij * (testDatNodes[rep(i, nrow(trainDatNodes)),] == trainDatNodes))
    w1 <- sapply(Time, function(x) sum(wij * (Time >= x)))
    exp(-sapply(Time, function(x)
      sum(Status * wij * (Time <= x) / w1, na.rm = TRUE)))
  })
  apply(predSV, 2, function(x) stepfun(sort(Time), c(1, x[order(Time)])))
}

#' Function used to compute the concordance measure at time tt
#'
#' @param tt is the time to evaluate the concordance measure
#' @param si is a list of survival functions output from getSurv()
#' @param sc is a Kaplan-Meier estimate of survival distribution of the censoring time
#' @param dat is a data frame for testing data
#'
#' This return a single numerical value
getCON <- function(tt, si, sc, dat) {
  c1 <- outer(dat$status, rep(1, length(dat$status)))
  c2 <- outer(dat$Time <= tt, dat$Time > tt)
  l <- length(si)
  stmp <- sapply(1:l, function(x) si[[x]](tt))
  c3 <- 1 * outer(stmp, stmp, "<") + 0.5 * outer(stmp, stmp, "==")
  c4 <- outer(sc(dat$Time + dat$D), sc(dat$D + tt))
  c4 <- ifelse(c4 < 1e-5, NA, c4)
  sum(c1 * c2 * c3 / c4, na.rm = TRUE) / sum(c1 * c2 / c4, na.rm = TRUE)
}

#' This return a single numerical value
getCON2 <- function(tt, si, sc, dat) {
  c1 <- outer(dat$status, rep(1, length(dat$status)))
  c2 <- outer(dat$Time <= tt, dat$Time > tt)
  l <- length(si)
  stmp <- sapply(1:l, function(x) si[[x]](tt))
  c3 <- 1 * outer(stmp, stmp, "<") + 0.5 * outer(stmp, stmp, "==")
  c4 <- outer(sc(dat$Time + dat$D), sc(dat$D + tt))
  c4 <- ifelse(c4 < 1e-5, NA, c4)
  n <- nrow(dat)
  sum(c1 * c2 * c3 / c4, na.rm = TRUE) / n / (n - 1)
}

#' Function used to compute the predicted survival without the ith subject using ensemble method
#'
#' @param fit is a ranger object
#' @param dat is a data frame for the training data (the data used in creating fit)
getSi <- function(fit, dat) {
  tt <- sort(dat$Time)
  rownames(dat) <- NULL
  sapply(1:nrow(dat), function(s) {
    fit.tmp <- fit
    takei <- which(sapply(fit$inbag.counts, function(x) x[s] == 0))
    fit.tmp$forest$num.trees <- fit.tmp$num.trees <- length(takei)
    fit.tmp$forest$child.nodeIDs <- fit.tmp$forest$child.nodeIDs[takei]
    fit.tmp$forest$split.values <- fit.tmp$forest$split.values[takei]
    fit.tmp$forest$chf <- fit.tmp$forest$chf[takei]
    fit.tmp$forest$split.varIDs <- fit.tmp$forest$split.varIDs[takei]
    fit.tmp$inbag.counts <- fit.tmp$inbag.counts[takei]
    getSurv(fit.tmp, dat, dat[s,])
  })
}
