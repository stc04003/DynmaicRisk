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
