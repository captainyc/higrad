library(higrad)

# get data from libsvm website; convert the data into data frame
rawdata <- read.table(url('https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/a9a'), 
                      sep = '\n', colClasses = 'character')
data <- data.frame(matrix(0, nrow(rawdata), 124))
for (i in 1:nrow(rawdata)) {
  row <- strsplit(rawdata[i, ], ' ')[[1]]
  data[i, 124] = as.integer(row[1])
  for (j in 2:length(row)) {
    temp <- as.integer(strsplit(row[j], ':')[[1]])
    data[i, temp[1]] = temp[2]
  }
}
colnames(data)[124] <- 'y'
# to save the data set for future use
# write.csv(data, file = 'adult_data.csv', row.names = FALSE)

# create the feature matrix x, and the response vector y
x <- data[, -ncol(data)]
y <- data$y
x <- scale(x)
x <- cbind(rep(1, nrow(x)), x)
colnames(x)[1] <- "X0"

# create training and test split
set.seed(1)
n_test <- 1000
index <- sample(1:nrow(x), n_test)
xtrain <- as.matrix(x[-index, ])
xtest <- as.matrix(x[index, ])
ytrain <- y[-index]
ytest <- y[index]

# parameters for higrad 
N <- nrow(xtrain) * 25
B <- 500
burnin <- N/5
n0 <- NA
eta <- 0.5
alpha <- 0.505
configs <- list()
# higrad
configs[[length(configs) + 1]] <- list(nsplits = 2, nthreads = 2, step.ratio = 1, n0 = n0, eta = eta, burnin = burnin, alpha = alpha)
# sgd
configs[[length(configs) + 1]] <- list(nsplits = 1, nthreads = 1, step.ratio = 1, n0 = 0, eta = eta, burnin = burnin, alpha = alpha)

# record for results
upper <- list()
lower <- list()
pred <- list()
beta <- list()
for (j in 1:length(configs)) {
  upper[[j]] <- matrix(0, B, n_test)
  lower[[j]] <- matrix(0, B, n_test)
  pred[[j]] <- matrix(0, B, n_test)
  beta[[j]] <- matrix(0, B, ncol(x))
}

# start experiment
for (b in 1:B) {
  for (j in 1:length(configs)) {
    fit <- higrad(xtrain, ytrain, model = "logistic", 
                  nsteps = N,
                  nsplits = configs[[j]]$nsplits, 
                  nthreads = configs[[j]]$nthreads, 
                  step.ratio = configs[[j]]$step.ratio,
                  n0 = configs[[j]]$n0,
                  burnin = configs[[j]]$burnin,
                  eta = configs[[j]]$eta, 
                  alpha = configs[[j]]$alpha,
                  replace = TRUE)
    prediction <- predict(fit, newx = xtest, alpha = 0.1, prediction.interval = TRUE)
    beta[[j]][b, ] <- coef(fit)
    pred[[j]][b, ] <- prediction$pred
    upper[[j]][b, ] <- prediction$upper
    lower[[j]][b, ] <- prediction$lower
  }
  if (b %% 10 == 0) {
    cat(b)
    save(pred, upper, lower, beta, file = "./results_adult.RData")
  }
}
save(pred, upper, lower, beta, file = "./results_adult.RData")

# plot
if (FALSE) {
  library(ggplot2)
  library(gridExtra)
  require(tikzDevice)
  options(tikzLatex = "/Library/TeX/texbin/pdflatex")
  
  # Figure 1, PI length vs average predicted probs
  j <- length(pred)
  ps <- 1/(1+exp(-pred[[j]]))
  p <- apply(ps, 2, median)
  u <- apply(ps, 2, function(x) quantile(x, 0.975))
  l <- apply(ps, 2, function(x) quantile(x, 0.025))
  pval <- apply(ps, 2, function(x) min(mean(x > 0.5), mean(x < 0.5))) * 2
  s <- u - l
  ind <- order(s, decreasing = TRUE)
  m <- 1000
  data.ci <- data.frame(id = 1:m, prob = p[ind[1:m]], upper = u[ind[1:m]], lower = l[ind[1:m]], pval = pval[ind[1:m]])
  
  tikz("./adult_data.tex", width = 4.5, height = 3)
  g3 <- ggplot(data.ci) +
    geom_point(aes(x = prob, y = upper - lower), size = 0.4) + 
    xlab("Average predicted probability") + 
    ylab("Empirical prediction interval length") +
    scale_y_log10(limits = c(1e-4, 1), breaks = c(0.0001, 0.001, 0.01, 0.1, 1), labels = c('0.01\\%', '0.1\\%', '1\\%', '10\\%', '100\\%')) +
    scale_x_continuous(labels = c('0\\%', '25\\%', '50\\%', '75\\%', '100\\%')) +
    theme_bw() + 
    theme(axis.text = element_text(size = 8), axis.title = element_text(size = 9))
  print(g3) 
  dev.off()
  
  # Figure 7, coverage plot
  B <- 100
  
  cover.higrad <- array(0, dim = c(B, B, 1000))
  for (b1 in 1:B) {
    for (b2 in 1:B) {
      cover.higrad[b1, b2, ] <- (pred[[1]][b1, ] <= upper[[1]][b2, ]) & (pred[[1]][b1, ] >= lower[[1]][b2, ])
    }
  }
  coverprob.higrad = apply(cover.higrad, 3, mean)
  
  cover.higrad <- array(0, dim = c(B, B, 1000))
  for (b1 in 1:B) {
    for (b2 in 1:B) {
      cover.higrad[b1, b2, ] <- (pred[[2]][b1, ] <= upper[[1]][b2, ]) & (pred[[2]][b1, ] >= lower[[1]][b2, ])
    }
  }
  coverprob.oracle = apply(cover.higrad, 3, mean)
  
  g1 <- ggplot(data.frame(coverage = coverprob.higrad), 
               aes(x = coverage)) + geom_histogram(binwidth = 0.05, color = "darkgray", fill = "gray") +
    scale_x_continuous(name = "Average coverage probability", breaks = seq(0, 1, 0.1), limits=c(0, 1)) +
    scale_y_continuous(name = "Count", limits = c(0, 450)) + 
    theme_bw() + 
    theme(axis.text=element_text(size = 6), axis.title=element_text(size = 8))
  
  g2 <- ggplot(data.frame(coverage = coverprob.oracle), 
               aes(x = coverage)) + geom_histogram(binwidth = 0.05, color = "darkgray", fill = "gray") +
    scale_x_continuous(name = "Average coverage probability", breaks = seq(0, 1, 0.1), limits=c(0, 1)) +
    scale_y_continuous(name = "Count", limits = c(0, 450)) +
    theme_bw() +
    theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8))
  
  tikz("./adult_coverage.tex", width = 6.2, height = 2.5)
  grid.arrange(g1, g2, ncol = 2)
  dev.off()
}