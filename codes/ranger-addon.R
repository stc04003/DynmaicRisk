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
  wbij <- sapply(1:fit$forest$num.trees, function(i) {
    out <- tabulate(1 + fit$forest$SampleIDs[[i]], nbins = n)
    return(out) })
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
#' @param df0 is a data frame for testing data
#'
#' This return a single numerical value
getCON <- function(tt, si, sc, df0) {
  c1 <- outer(df0$status, rep(1, length(df0$status)))
  c2 <- outer(df0$Time <= tt, df0$Time > tt)
  l <- length(si)
  stmp <- sapply(1:l, function(x) si[[x]](tt))
  c3 <- 1 * outer(stmp, stmp, "<") + 0.5 * outer(stmp, stmp, "==")
  c4 <- outer(sc(df0$Time + df0$D), sc(df0$Time + tt))
  c4 <- ifelse(c4 < 1e-5, NA, c4)
  sum(c1 * c2 * c3 / c4, na.rm = TRUE) / sum(c1 * c2 / c4, na.rm = TRUE)
}


#' Function used to compute the predicted survival without the ith subject using ensemble method
#'
#' @param fit is a ranger object
#' @param df is a data frame for the training data (the data used in creating fit)
getSi <- function(fit, df) {
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
