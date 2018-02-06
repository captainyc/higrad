# Simulations to compare HiGrad's estimation accuracy with multiple configurations and settings.
# We also compare the performance to a vanilla SGD

library(higrad)
library(optparse)

# parse command line argument
option_list = list(
  make_option(c('-m', '--model'), action = 'store', default = 'lm', type = 'character',
              help = 'regression type, linear or logistic'),
  make_option(c('-d', '--dimension'), action = 'store', default = 50, type = 'numeric',
              help = 'dimension of the parameter'),
  make_option(c('-r', '--num.repeat'), action = 'store', default = 100, type = 'numeric',
              help = 'number of repeats'),
  make_option(c('-t', '--theta'), action = 'store', default = 'dense', type = 'character',
              help = 'type of true theta')
)
opt <- parse_args(OptionParser(option_list=option_list))

# set up parameters
model <- opt$model
d <- opt$dimension
num.repeat <- opt$num.repeat
if (opt$theta == 'uniform') {
  theta.star <- seq(0, 1, length = d)
} else if (opt$theta == 'null') {
  theta.star <- rep(0, d)
} else if (opt$theta == 'dense') {
  theta.star <- rep(1/sqrt(d), d)
} else if (opt$theta == 'sparse') {
  theta.star <- c(rep(1, ceiling(d/10)), rep(0, d - ceiling(d/10)))
  theta.star <- theta.star / sqrt(sum(theta.star^2))
} else {
  print('theta type not supported. Using null hypothesis...')
  theta.star <- rep(0, d) 
}
Ns <- round(10^seq(4, 6, length = 11))
sigma <- 1
alpha <- 0.55
eta <- ifelse(model == 'lm', 0.1, 0.4)
burnin <- 0

filename <- paste(model, 'd', d, 'theta', opt$theta, sep = '_')

# print simulation type
cat('Regression type:', model, '\n')
cat('Dimension:', d, '\n')
cat('Starting simulation...\n')

# create configurations
# 4 configs considered (4, 1-4, 2-4, 1-2-4)
configs <- list()
for (n0 in c(0, NA)) {
  configs[[length(configs) + 1]] <- list(nsplits = 1, nthreads = 4, step.ratio = 1, n0 = n0, eta = eta, burnin = burnin)
  configs[[length(configs) + 1]] <- list(nsplits = 2, nthreads = 2, step.ratio = 1, n0 = n0, eta = eta, burnin = burnin)
}
configs[[length(configs) + 1]] <- list(nsplits = 1, nthreads = 1, step.ratio = 1, n0 = 0, eta = eta, burnin = burnin)

record <- list(length(configs))

# create record list
for (j in 1:length(configs)) {
  record[[j]] <- list(estimate = array(0, dim = c(num.repeat, d, length(Ns))))
}

# start simulation
for (k in 1:length(Ns)) {
  N <- Ns[k]
  for (i in 1:num.repeat) {
    set.seed(i)
    x <- matrix(rnorm(N * d), N, d)
    if (model == 'lm') {
      y <- as.numeric(x %*% theta.star) + rnorm(N, 0, sigma)
    } else {
      y <- ifelse(runif(N, 0, 1) > 1 / (1 + exp(-as.numeric(x %*% theta.star))), -1, 1)
    }
    
    for (j in 1:length(configs)) {
      fit <- higrad(x, y, model = model, 
                    nsplits = configs[[j]]$nsplits, 
                    nthreads = configs[[j]]$nthreads, 
                    step.ratio = configs[[j]]$step.ratio,
                    n0 = configs[[j]]$n0,
                    burnin = configs[[j]]$burnin,
                    eta = configs[[j]]$eta, 
                    alpha = alpha)
      record[[j]]$estimate[i, , k] <- fit$coefficients
    }
  }
}

save(record, file = paste0('./record_accuracy_', filename, '.RData'))