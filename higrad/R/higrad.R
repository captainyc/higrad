#' @importFrom grDevices rainbow
#' @importFrom graphics axis
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom stats qt
#' @importFrom stats rnorm


#' @title Fitting HiGrad
#'
#' @description
#' \code{higrad} is used to implement hierarchical incremental gradient descent (HiGrad), an algorithm that conducts statistical inference for online learning.
#'
#' @details
#' HiGrad is designed to conduct statistical inference for online learning, without incurring additional computational cost compared with the vanilla stochastic gradient descent (SGD).
#' The HiGrad procedure begins by performing SGD iterations for a while and then split the single thread into a few, and this procedure hierarchically operates in this fashion along each thread.
#' With predictions provided by multiple threads in place, a t-based confidence interval is constructed by de-correlating predictions using covariance structures given by the Ruppert–Polyak averaging scheme.
#' In order to implement HiGrad, a configuration of the tree structure needs to be specified. The default setting is a binary tree with 2 splits.
#' The step size is set to be \code{eta*t^(-alpha)}.
#'
#' @param x input matrix of features. Each row is an observation vector, and each column is a feature.
#' @param y response variable. Quantitative for \code{model = "lm"}. For \code{model = "logistic"} it should be a factor with two levels.
#' @param model type of model to fit. Currently only linear regression (\code{"lm"}) and logistic regression (\code{"logistic"}) are supported.
#' @param nsteps total number of steps. This is equivalent to the number of queries made to get a noisy evaluation of the gradient.
#' @param nsplits number of splits in the HiGrad tree.
#' @param nthreads numbers of threads each previous thread is split into. Either a number (equal split size throughout) or a vector.
#' @param step.ratio ratio of the lengths of the threads from the two adjacent levels (the latter one divided by the previous). Either a number (equal ratio throughout) or a vector.
#' @param n0 length of the 0th-level thread.
#' @param skip number of steps to skip when estimating the coefficients by averaging.
#' @param eta constant in front of the step size. See Details for the formula of the step size.
#' @param alpha exponent of the step size. See Details for the formula of the step size.
#' @param burnin number of steps as the burn-in period. The burn-in period is not accounted for in the total budget \code{nsteps}.
#' @param start starting values of the coefficients.
#' @param replace logical; whether or not to sample the data with replacement.
#' @param track logical; whether or not to store the entire path for plotting.
#'
#' @references Weijie Su and Yuancheng Zhu. (2018) \emph{Statistical Inference for Online Learning and Stochastic Approximation via Hierarchical Incremental Gradient Descent}. \url{https://arxiv.org/abs/1802.04876}.
#' @return An object with S3 class \code{higrad}.
#' \item{coefficients}{estimate of the coefficients.}
#' \item{coefficients.bootstrap}{matrix of estimates of the coefficients along each HiGrad threads.}
#' \item{model}{model type.}
#' \item{Sigma0}{covariance structure \eqn{\Sigma} of the estimates.}
#' \item{track}{entire path of the estimates along each thread. Can be used for diagnostic and check for convergence.}
#'
#' @seealso See \code{\link{print.higrad}}, \code{\link{plot.higrad}}, \code{\link{predict.higrad}} for other methods for the \code{higrad} class.
#'
#' @examples
#' # fitting linear regression on a simulated dataset
#' n <- 1e3
#' d <- 10
#' sigma <- 0.1
#' theta <- rep(1, d)
#' x <- matrix(rnorm(n * d), n, d)
#' y <- as.numeric(x %*% theta + rnorm(n, 0, sigma))
#' fit <- higrad(x, y, model = "lm")
#' print(fit)
#' # predict for 10 new samples
#' newx <- matrix(rnorm(10 * d), 10, d)
#' pred <- predict(fit, newx)
#' pred
#'
#' @export
higrad <- function(x, y, model = "lm",
                   nsteps = nrow(x), nsplits = 2, nthreads = 2, step.ratio = 1, n0 = NA, skip = 0,
                   eta = 1/2, alpha = 1/2, burnin = round(nsteps / 10), start = rnorm(ncol(x), 0, 0.01),
                   replace = FALSE, track = FALSE) {
  # check for missing or unsupported inputs
  if (missing(x)) {
    stop("'x' not specified")
  }
  if (missing(y)) {
    stop("'y' not specified")
  }
  if (!(model %in% c("lm", "logistic"))) {
    stop("'model' not supported")
  }
  if (model == "logistic") {
    if (length(unique(y)) != 2) stop("response is not a binary variable")
    else y <- as.numeric(factor(y)) * 2 - 3
  }
  x <- as.matrix(x)
  if (nrow(x) < 2) {
    stop("'x' should be a matrix with 2 or more columns")
  }
  if (nrow(x) != length(y)) {
    stop("the dimensions of 'x' and 'y' do not agree")
  }

  # create split parameters
  split <- createSplit(nsteps, nsplits, nthreads, step.ratio, n0)
  ns <- split$ns
  K <- split$K
  Bs <- split$Bs
  n0 <- split$n0
  B <- prod(Bs)
  ws <- c(n0, ns) * cumprod(c(1, Bs))
  ws <- ws / sum(ws)
  d <- ncol(x)

  # create gradient function
  if (model == "lm") {
    getGradient <- function(theta, x1, y1) { x1 * (sum(theta * x1) - y1) }
  } else {
    getGradient <- function(theta, x1, y1) { -y1 * x1 / (1 + exp(y1 * sum(theta * x1))) }
  }

  # create step function
  stepSize <- function(n) { eta / n^alpha }

  # create sample function
  sampleNext <- function(current, n, replace) { ifelse(replace, sample(n, 1), current %% n + 1) }
  idx <- 0

  # fit higrad
  # burnin stage
  if (burnin > 0) {
    for (i in 1:burnin) {
      idx <- sampleNext(idx, nrow(x), replace)
      start <- start - stepSize(i) * getGradient(start, x[idx, ], y[idx])
    }
  }

  # create a matrix that stores stagewise average
  theta.avg <- array(NA, dim = c(B, d, K+1))
  # set up theta.track for plotting
  if (track) theta.track <- matrix(start, d, 1)

  # zeroth stage
  # theta matrix for the current stage
  theta.current <- matrix(NA, d, n0+1)
  # initial value
  theta.current[, 1] <- start
  # iteration
  if (n0 > 0) {
    for (i in 1:n0) {
      idx <- sampleNext(idx, nrow(x), replace)
      theta.current[, i+1] <- theta.current[, i] - stepSize(i) * getGradient(theta.current[, i], x[idx, ], y[idx])
    }
    # average and store in theta.avg
    theta.avg[, , 1] <- matrix(rowMeans(theta.current[, -(1:(floor(n0 * skip)+1)), drop = FALSE]), B, d, byrow = TRUE)
  } else {
    theta.avg[, , 1] <- matrix(start, B, d, byrow = TRUE)
  }
  # concatenate theta.track
  if (track) theta.track <- cbind(theta.track, theta.current[, -1])
  # set initial value for the next stage
  start <- matrix(rep(theta.current[, n0+1], each = Bs[1]), Bs[1], d)

  # iterative through stages
  for (k in 1:K) {
    # theta matrix for the current stage
    theta.current <- array(NA, dim = c(cumprod(Bs)[k], d, ns[k]+1))
    # initial value
    theta.current[, , 1] <- start
    # iteration
    n.current <- cumsum(c(n0, ns))[k]
    for (i in 1:ns[k]) {
      for (j in 1:cumprod(Bs)[k]) {
        idx <- sampleNext(idx, nrow(x), replace)
        theta.current[j, , i+1] <- theta.current[j, , i] - stepSize(n.current+i) * getGradient(theta.current[j, , i], x[idx, ], y[idx])
      }
    }
    # average and store in theta.avg
    theta.avg[, , k+1] <- matrix(rep(apply(theta.current[, , -(1:(floor(ns[k] * skip)+1)), drop = FALSE], 1:2, mean), each = B / cumprod(Bs)[k]), B, d)
    # concatenate theta.track
    if (track) theta.track <- cbind(theta.track, theta.current[1, , -1])
    # set initial value for the next stage
    if (k < K) {
      start <- matrix(rep(theta.current[, , ns[k]+1], each = Bs[k+1]), cumprod(Bs)[k+1], d)
    }
  }

  # weighted average across stages
  thetas <- t(sapply(1:B, function(i) theta.avg[i, , ] %*% ws))

  out <- list()
  out$coefficients <- colMeans(thetas)
  out$coefficients.bootstrap <- thetas
  out$model <- model
  out$Sigma0 <- getSigma0(c(n0, ns), c(1, Bs), ws)
  out$track <- ifelse(track, theta.track, NA)

  class(out) <- "higrad"

  out
}

#' @title Obtain Prediction and Confidence Intervals From a HiGrad Fit
#'
#' @description \code{predict} can be applied with a \code{higrad} object to obtain predictions and confidence intervals.
#'
#' @param object a fitted object of class \code{higrad}.
#' @param newx matrix of new values for \code{x} at which predictions are to be made.
#' @param alpha significance level. The confidence level of the interval is thus 1 - \code{alpha}.
#' @param type type of prediction required. Type \code{"link"} gives the linear predictors for \code{"logistic"}; for \code{"lm"} models it gives the fitted values. Type \code{"response"} gives the fitted probabilities for \code{"logistic"}; for \code{"lm"} type \code{"response"} is equivalent to type \code{"link"}.
#' @param prediction.interval logical; indicator of whether prediction intervals should be returned instead of confidence intervals.
#' @param ... other prediction options.
#'
#' @return A list with components
#' \item{pred}{predicted values.}
#' \item{upper}{upper limit of the confidence/prediction intervals.}
#' \item{lower}{lower limit of the confidence/prediction intervals.}
#'
#' @export
predict.higrad <- function(object, newx, alpha = 0.05, type = "link", prediction.interval = FALSE, ...) {
  if (is.vector(newx)) {
    newx = matrix(newx, 1, length(newx))
  }
  if (ncol(newx) != length(object$coefficients)) {
    stop("'newx' has the wrong dimension")
  }
  Sigma0 <- object$Sigma0
  B <- nrow(Sigma0)
  mu <- as.numeric(newx %*% object$coefficients)

  # n x B matrix of predicted values
  mus <- newx %*% t(object$coefficients.bootstrap)
  # standard errors
  ses <- sqrt(sum(Sigma0) * colSums(t(mus - mu) * solve(Sigma0, t(mus - mu))) / (B^2 * (B-1)))
  if (prediction.interval) {
    margin <- qt(1-alpha/2, B-1) * ses * sqrt(2)
  } else {
    margin <- qt(1-alpha/2, B-1) * ses
  }
  upper <- mu + margin
  lower <- mu - margin

  if (object$model == "logistic" & type == "response") {
    mu = 1 / (1 + exp(-mu))
    upper = 1 / (1 + exp(-upper))
    lower = 1 / (1 + exp(-lower))
  }

  out <- list(pred = mu, upper = upper, lower = lower)
  return(out)
}

#' @title Print a \code{higrad} Object
#'
#' @description Print the coefficients estimates obtained by HiGrad.
#'
#' @param x a fitted object of class \code{higrad}.
#' @param ... additional print arguments.
#'
#' @export
print.higrad <- function(x, ...) {
  print(x$coefficients, ...)
}

#' @title Plot a \code{higrad} Object
#'
#' @description Produces a coefficient paths for a fitted \code{higrad} object.
#'
#' @param x a fitted object of class \code{higrad}.
#' @param ... additional graphical parameters.
#'
#' @export
plot.higrad <- function(x, ...) {
  if (!is.na(x$track)) {
    track <- x$track
    n <- ncol(track)
    d <- nrow(track)
    colors <- rainbow(d)
    n.subsample <- round(seq(1, n, length = 100))
    track <- track[, n.subsample]
    plot(n.subsample, track[1, ], type = "n", ylim = c(min(track), max(track)), xlab = "Step", ylab = "Coefficient estimates")
    for (i in 1:d) {
      lines(n.subsample, track[i, ], col = colors[i])
    }
    axis(side = 1)
    axis(side = 2)
  }
}

# create configurations based on the input parameters
createSplit <- function(nsteps, nsplits, nthreads, step.ratio, n0) {
  N <- nsteps
  K <- nsplits
  if (length(nthreads) == 1) {
    Bs <- rep(nthreads, K)
  } else {
    Bs <- nthreads
  }
  stopifnot(length(Bs) == K)
  step.ratio <- step.ratio^(1:K)
  if (is.na(n0)) {
    ns <- round((N / sum(c(1, step.ratio) * c(1, cumprod(Bs)))) * c(1, step.ratio))
    n0 <- ns[1]
    ns <- ns[-1]
  } else {
    ns <- round(((N - n0) / sum(step.ratio * cumprod(Bs))) * step.ratio)
  }
  ns <- ns[ns > 0]
  K <- sum(ns > 0)
  Bs <- Bs[ns > 0]

  return(list(ns = ns, K = K, Bs = Bs, n0 = n0))
}

# generate the covariance matrix of the estimates along different threads
getSigma0 <- function(ns, Bs, ws = NULL) {
  B <- prod(Bs)
  BBs <- cumprod(Bs)
  if (is.null(ws)) {
    ws <- ns * BBs
    ws <- ws / sum(ws)
  }
  Sigma0 <- matrix(0, B, B)
  for (l in which(ns > 0)[1]:length(ns)) {
    Sigma0 <- Sigma0 + Matrix::bdiag(rep(list(matrix(ws[l]^2/ns[l], B/BBs[l], B/BBs[l])), BBs[l]))
  }
  return(as.matrix(Sigma0))
}
