
## Dynamic risk prediction triggered by intermediate events using survival tree ensembles

This repository contains example codes used to run the simulation in the
manuscript “Dynamic risk prediction triggered by intermediate events
using survival tree ensembles” described in
<https://arxiv.org/abs/2011.07175>.

## Example

In the following example, we used the `ranger()` function implemented in
`R` package `ranger` (Wright and Ziegler, 2017) to construct survival
trees then apply the proposed ensemble method to compute the predicted
survival probabilities.

We consider a scenario where the event times are generated from an
irreversible multi-state model with three states: healthy, diseased, and
death. We assume that all subjects started in the healthy state, disease
onset is an intermediate event, and death is the event of interest.
Suppose we are interested in predicting the survival probability among
those who experienced an intermediate event, the target landmark
probability is
![equation](https://latex.codecogs.com/svg.image?P\(T%5Cge%20U%20+%20t%7CT%5Cge%20U,%20U,%20W\(U\),%20Z\),%20)
where ![equation](https://latex.codecogs.com/svg.image?T) is a
continuous failure time,
![equation](https://latex.codecogs.com/svg.image?U) is the time to the
intermediate event (time from the healthy state to the disease state),
![equation](https://latex.codecogs.com/svg.image?W\(U\)) is a
q-dimensional vector of time-dependent predictors measured at
![equation](https://latex.codecogs.com/svg.image?U), and
![equation](https://latex.codecogs.com/svg.image?Z) is a p-dimensional
vector of baseline predictors.

The following codes generate simulated landmark data from the
irreversible multi-state model with the intermediate event as the
landmark time.

``` r
> source("codes/data.R")
> set.seed(1); dat <- simDat(n = 400, cen = .2)
```

The `simDat()` function from `codes/data.R` takes on two arguments, `n`
and `cen` for sample size and censoring percentage at the baseline,
respectively.

``` r
> names(dat)
```

``` 
 [1] "Time"   "status" "D"      "Z.1"    "Z.2"    "Z.3"    "Z.4"    "Z.5"   
 [9] "Z.6"    "Z.7"    "Z.8"    "Z.9"    "Z.10"   "W.1"    "W.2"    "W.3"   
[17] "W.4"    "W.5"    "W.6"    "W.7"    "W.8"    "W.9"    "W.10"  
```

``` r
> head(dat[, c(1:5, 14:15)]) 
```

``` 
        Time status         D        Z.1        Z.2           W.1           W.2
1 0.33287199      1 3.2491162 -0.6264538  1.0744410 -0.2560830398 -0.2965939716
2 8.41871360      1 1.0447505  0.1836433  1.8956548 -0.0318106031 -0.0070139599
3 0.13400898      1 4.4229831 -0.8356286 -0.6029973  0.3014749450  0.3778652986
4 1.09821920      1 2.3055208  1.5952808 -0.3908678 -0.0305878910  0.1440122438
5 0.02807212      1 0.1503251  0.3295078 -0.4162220  0.0003901334 -0.0008796919
8 0.58646119      0 2.9984051  0.7383247 -0.2956775 -0.1781971479  0.1950746982
```

The `simDat()` returns a `data.frame` with the following variables.

  - `Time` observed survival time after
    ![equation](https://latex.codecogs.com/svg.image?U), e.g.,
    ![equation](https://latex.codecogs.com/svg.image?T&space;-&space;U).
  - `status` censoring indicator; 1 if `Time` is the event of interest
    (death) and 0 if `Time` is censored.
  - `D` time from the healthy state to the diseased state. This is also
    the landmark time in this scenario.
  - `Z.1`…`Z.10` are 10 baseline predictors.
  - `W.1`…`W.10` are 10 time-dependent predictors measured at
    ![equation](https://latex.codecogs.com/svg.image?U).

The following fits a random survival forest with `ranger()`.

``` r
> library(survival)
> library(ranger)
> fit <- ranger(Surv(Time, status) ~ ., data = dat)
```

The proposed ensemble procedure is implemented in `getSurv()` from
`codes/ranger-addon.R`. The following codes generate a testing data and
computes the predicted probabilities with the proposed ensemble
procedure.

``` r
> source("codes/ranger-addon.R")
> dat0 <- simDat(500, 0)
> predS <- getSurv(fit, dat, dat0)
> head(predS, 3)
```

    [[1]]
    Step function
    Call: stepfun(sort(Time), c(1, x[order(Time)]))
     x[1:239] = 0.0081491, 0.011054, 0.011668,  ..., 76.389, 113.75
    240 plateau levels =      1,      1, 0.99952,  ..., 0.091544, 0.091544
    
    [[2]]
    Step function
    Call: stepfun(sort(Time), c(1, x[order(Time)]))
     x[1:239] = 0.0081491, 0.011054, 0.011668,  ..., 76.389, 113.75
    240 plateau levels =      1,  0.999,  0.999,  ..., 0.028241, 0.028241
    
    [[3]]
    Step function
    Call: stepfun(sort(Time), c(1, x[order(Time)]))
     x[1:239] = 0.0081491, 0.011054, 0.011668,  ..., 76.389, 113.75
    240 plateau levels =      1, 0.99801, 0.99801,  ..., 0.046578, 0.046578

The `getSurv()` function returns a list, where the ith member of the
list gives the landmark probability for the ith subject from the testing
data. Due to the complicated relationship between event times and
longitudinal markers, deriving the closed-form expression of the true
probability is challenging. The `getTrueSurv()` function from
`codes/data.R` compute the true landmark probability using Monte Carlo
method for this scenario.

``` r
> trueS <- getTrueSurv(dat0)
```

The `getTrueSurv()` function returns a list, where the ith member of the
list gives the true landmark probability for the ith subject from the
testing data. The following codes calculate the integrated mean absolute
error and the integrated mean squared error to evaluate the model
performance.

``` r
> t0 <- seq(0, quantile(dat0$Time, .9), .01)
> ## Intergrated absolute error
> mean(sapply(1:nrow(dat0), function(i) abs(predS[[i]](t0) - trueS[[i]](t0))))
```

    [1] 0.1817185

``` r
> ## Intergrated MSE
> mean(sapply(1:nrow(dat0), function(i) (predS[[i]](t0) - trueS[[i]](t0))^2))
```

    [1] 0.06385354

The considered example is the Scenario (III-A) in our manuscript. All
codes are available in the `code` folder.

## Reference

  - Wright, M. N. and Ziegler, A. (2017). `ranger`: A fast
    implementation of random forests for high dimensional data in `C++`
    and `R`. Journal of Statistical Software 77 1–17.
