<!-- README.md is generated from README.Rmd. Please edit that file -->
higrad
======

The goal of higrad is to implement the Hierarchical Incremental GRAdient Descent (HiGrad) algorithm. HiGrad is a first-order algorithm for finding the minimizer of a function in online learning just like SGD and, in addition, this method attaches a confidence interval to assess the uncertainty of its predictions.

Example
-------

This is a basic example which shows you how to solve a linear regression using higrad with simulated data. The predictions obtained at the end come with 95% confidence intervals.

``` r
library(higrad)
# generate a data set for linear regression
n <- 1e6
d <- 50
sigma <- 1
theta <- rep(0, d)
x <- matrix(rnorm(n * d), n, d)
y <- as.numeric(x %*% theta + rnorm(n, 0, sigma))
# fit the linear regression with higrad using the default setting
fit <- higrad(x, y, model = "lm")
# predict for 10 new samples
newx <- matrix(rnorm(10 * d), 10, d)
pred <- predict(fit, newx)
```
