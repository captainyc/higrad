# R code for reproducing the plots in Su & Zhu (2018)

library(higrad)
require(xtable)
require(ggplot2)
require(grid)
require(gridExtra)
require(tikzDevice)
options(tikzLatex="/Library/TeX/texbin/pdflatex")

###-----------------------------------------------------------------
# Figure 5

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

Ns <- round(10^seq(4, 6, length = 11))
B <- 100
d <- 50

models <- c("lm", "logistic")
theta.types <- c("null", "sparse", "dense")

tikz("./figs/fig_accuracy.tex", height = 6, width = 11)
g <- list()
for (model in models) {
  for (theta.type in theta.types) {
    if (theta.type == "null") {
      theta.star <- rep(0, d)
    } else if (theta.type == "dense") {
      theta.star <- rep(1/sqrt(d), d)
    } else if (theta.type == "sparse") {
      theta.star <- c(rep(1, ceiling(d/10)), rep(0, d - ceiling(d/10)))
      theta.star <- theta.star / sqrt(sum(theta.star^2))
    } else {
      theta.star <- seq(0, 1, length = d)
    }
    
    filename <- paste(model, 'd', d, 'theta', theta.type, sep = '_')
    load(paste0("./record_accuracy_", filename, '.RData'))
    risk <- matrix(0, length(record), length(Ns))
    for (j in 1:length(record)) {
      risk[j, ] <- sapply(1:length(Ns), function(k) mean(sqrt(rowMeans((record[[j]]$estimate[1:B, , k] - matrix(theta.star, B, d, byrow = TRUE))^2))))
    }
    r <- risk / matrix(risk[length(record),], length(record), length(Ns), byrow = TRUE)
    temp <- data.frame(N = rep(Ns, 4), risk = c(r[1, ], r[3, ], r[5, ], r[7, ]), K = factor(rep(c(2, 4), each = length(Ns) * 2)), N0 = factor(rep(c(TRUE, FALSE, TRUE, FALSE), each = length(Ns))))
    g[[length(g)+1]] <- ggplot(temp, aes(x = N, y = risk, colour = interaction(K, N0), linetype = interaction(K, N0), group = interaction(K, N0))) + 
      scale_color_manual(values=c("black", "red", "black", "red")) + 
      scale_linetype_manual(values = c("dashed", "dotdash", "solid", "dotted")) + 
      scale_x_log10() + 
      scale_y_continuous(limits = c(1, 1.35)) + 
      ylab("Normalized risk") + 
      xlab("Total number of steps") +
      geom_line() + 
      theme_bw() +
      theme(legend.position="none")+
      ggtitle(paste(ifelse(model == "lm", "Linear regression", "Logistic regression"), theta.type, sep = ", "), subtitle = NULL)
  }
}
multiplot(g[[1]], g[[4]], g[[2]], g[[5]], g[[3]], g[[6]], cols=3)
dev.off()

# legend
tikz("./figs/legend_accuracy.tex", height = 0.25, width = 5)
par(mfrow = c(1, 4), mar = rep(0, 4))
l = 0.5
a = 1/4
b = 0.5
plot(0, 0, "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-l, 1/2), ylim = c(-1, 1))
segments(rep(0, 4), c(-1.5, -0.5, 0.5, 1.5)*b, rep(a, 4), c(-1.5, -0.5, 0.5, 1.5)*b)
segments(0, -1.5*b, 0, 1.5*b, lty = 3, col = "gray")
legend("left", lty = "solid", lwd = 1, legend = ":", bty = "n", seg.len = 4)
a = 1/5
b = 0.5
plot(0, 0, "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-l, 1/2), ylim = c(-1, 1))
segments(0, 0, a, 0)
segments(rep(a, 4), c(-1.5, -0.5, 0.5, 1.5)*b, rep(2*a, 4), c(-1.5, -0.5, 0.5, 1.5)*b)
segments(a, -1.5*b, a, 1.5*b, lty = 3, col = "gray")
legend("left", lty = "dotted", col = "red", lwd = 1, legend = ":", bty = "n", seg.len = 4)
a = 1/6
b = 0.5
plot(0, 0, "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-l, 1/2), ylim = c(-1, 1))
segments(rep(0, 2), c(-b, b), rep(a, 2), c(-b, b))
segments(rep(a, 4), c(-1.5, -0.5, 0.5, 1.5)*b, rep(2*a, 4), c(-1.5, -0.5, 0.5, 1.5)*b)
segments(0, -b, 0, b, lty = 3, col = "gray")
segments(a, -1.5*b, a, -0.5*b, lty = 3, col = "gray")
segments(a, 0.5*b, a, 1.5*b, lty = 3, col = "gray")
legend("left", lty = "dashed", lwd = 1, legend = ":", bty = "n", seg.len = 4)
a = 1/7
b = 0.5
plot(0, 0, "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-l, 1/2), ylim = c(-1, 1))
segments(0, 0, a, 0)
segments(rep(a, 2), c(-b, b), rep(2*a, 2), c(-b, b))
segments(rep(2*a, 4), c(-1.5, -0.5, 0.5, 1.5)*b, rep(3*a, 4), c(-1.5, -0.5, 0.5, 1.5)*b)
segments(a, -b, a, b, lty = 3, col = "gray")
segments(2*a, -1.5*b, 2*a, -0.5*b, lty = 3, col = "gray")
segments(2*a, 0.5*b, 2*a, 1.5*b, lty = 3, col = "gray")
legend("left", lty = "dotdash", col = "red", lwd = 1, legend = ":", bty = "n", seg.len = 4)
dev.off()

###-----------------------------------------------------------------
# Figure 6
models <- c("lm", "logistic")
theta.types <- c("null", "sparse", "dense")
d <- 50
N <- 1e6
B <- 100
n0s <- c(0, NA)

results.higrad <- data.frame(dim = integer(), 
                             N = integer(),
                             model = factor(levels = c("lm", "logistic")), 
                             theta.type = factor(levels = theta.types),
                             config = integer(), 
                             coverage = double(),
                             ci.length = double(),
                             estimate = double())

configs <- list()

for (n0 in c(0, NA)) {
  configs[[length(configs) + 1]] <- list(nsplits = 1, nthreads = 2, step.ratio = 1, n0 = n0)
  configs[[length(configs) + 1]] <- list(nsplits = 1, nthreads = 4, step.ratio = 1, n0 = n0)
  configs[[length(configs) + 1]] <- list(nsplits = 2, nthreads = 2, step.ratio = 1, n0 = n0)
}
configs[[length(configs) + 1]] <- list(nsplits = 1, nthreads = 1, step.ratio = 1, n0 = 0)

for (model in models) {
  for (theta.type in theta.types) {
    if (theta.type == "null") {
      theta.star <- rep(0, d)
    } else if (theta.type == "dense") {
      theta.star <- rep(1/sqrt(d), d)
    } else if (theta.type == "sparse") {
      theta.star <- c(rep(1, ceiling(d/10)), rep(0, d - ceiling(d/10)))
      theta.star <- theta.star / sqrt(sum(theta.star^2))
    } else {
      theta.star <- seq(0, 1, length = d)
    }
    
    filename <- paste(model, 'd', d, 'theta', theta.type, sep = '_')
    if (file.exists(paste0('./record_', filename, '.RData'))){
      load(paste0('./record_coverage_', filename, '.RData'))
    }
    for (j in 1:length(configs)) {
      results.higrad[nrow(results.higrad)+1, ] = data.frame(
        d, 
        N,
        model,
        theta.type,
        j,
        round(mean(record[[j]]$coverage[1:B, ]), 4),
        round(mean(record[[j]]$ci.length[1:B, ]), 4),
        round(mean(sqrt(rowSums((record[[j]]$estimate[1:B, ] - matrix(theta.star, B, d, byrow = TRUE))^2))), 4)
      )
    }
  }
}

MakePlot <- function(d, reg.type, tikz = TRUE) {
  filename <- paste0('./figs/fig_coverage_', reg.type, '.tex')
  
  data <- subset(results.higrad, dim == d & model == reg.type & config < 7)
  data <- data[order(data$config), ]
  data$number <- as.factor(nrow(data):1)
  
  gray = gray.colors(3, start = 0.5, end = 0.9)
  
  g.mid = ggplot() +
    scale_x_continuous(limits = c(0, 0.7)) +
    scale_y_continuous(limits = c(-3, 3), expand = c(0.036, 0.036)) +
    geom_segment(aes(x = 0, xend = l/2, y = h*2 - u*2.5, yend = h*2 - u*2.5)) + 
    geom_segment(aes(x = 0, xend = l/2, y = -h*2 - u*2.5, yend = -h*2 - u*2.5)) +
    geom_segment(aes(x = 0, xend = 0, y = -h*2 - u*2.5, yend = h*2 - u*2.5), linetype = "dashed", color = "gray") +
    ###
    geom_segment(aes(x = 0, xend = l/4, y = h - u*1.5, yend = h - u*1.5)) + 
    geom_segment(aes(x = 0, xend = l/4, y = -h - u*1.5, yend = -h - u*1.5)) +
    geom_segment(aes(x = 0, xend = l/4, y = h*3 - u*1.5, yend = h*3 - u*1.5)) + 
    geom_segment(aes(x = 0, xend = l/4, y = -h*3 - u*1.5, yend = -h*3 - u*1.5)) +
    geom_segment(aes(x = 0, xend = 0, y = -h*3 - u*1.5, yend = h*3 - u*1.5), linetype = "dashed", color = "gray") +
    ###
    geom_segment(aes(x = 0, xend = l/6, y = h*2 - u*.5, yend = h*2 - u*.5)) + 
    geom_segment(aes(x = 0, xend = l/6, y = -h*2 - u*.5, yend = -h*2 - u*.5)) +
    geom_segment(aes(x = l/6, xend = l/3, y = h - u*.5, yend = h - u*.5)) + 
    geom_segment(aes(x = l/6, xend = l/3, y = -h - u*.5, yend = -h - u*.5)) +
    geom_segment(aes(x = l/6, xend = l/3, y = h*3 - u*.5, yend = h*3 - u*.5)) + 
    geom_segment(aes(x = l/6, xend = l/3, y = -h*3 - u*.5, yend = -h*3 - u*.5)) +
    geom_segment(aes(x = 0, xend = 0, y = -h*2 - u*.5, yend = h*2 - u*.5), linetype = "dashed", color = "gray") +
    geom_segment(aes(x = l/6, xend = l/6, y = -h*3 - u*.5, yend = h*3 - u*.5), linetype = "dashed", color = "gray") +
    ###
    geom_segment(aes(x = 0, xend = l/3, y = u*0.5, yend = u*0.5)) + 
    geom_segment(aes(x = l/3, xend = l*2/3, y = h*2 + u*0.5, yend = h*2 + u*0.5)) + 
    geom_segment(aes(x = l/3, xend = l*2/3, y = -h*2 + u*0.5, yend = -h*2 + u*0.5)) +
    geom_segment(aes(x = l/3, xend = l/3, y = -h*2 + u*0.5, yend = h*2 + u*0.5), linetype = "dashed", color = "gray") +
    ###
    geom_segment(aes(x = 0, xend = l/5, y = u*1.5, yend = u*1.5)) + 
    geom_segment(aes(x = l/5, xend = l*2/5, y = h + u*1.5, yend = h + u*1.5)) + 
    geom_segment(aes(x = l/5, xend = l*2/5, y = -h + u*1.5, yend = -h + u*1.5)) +
    geom_segment(aes(x = l/5, xend = l*2/5, y = h*3 + u*1.5, yend = h*3 + u*1.5)) + 
    geom_segment(aes(x = l/5, xend = l*2/5, y = -h*3 + u*1.5, yend = -h*3 + u*1.5)) +
    geom_segment(aes(x = l/5, xend = l/5, y = -h*3 + u*1.5, yend = h*3 + u*1.5), linetype = "dashed", color = "gray") +
    ###
    geom_segment(aes(x = 0, xend = l/7, y = u*2.5, yend = u*2.5)) + 
    geom_segment(aes(x = l/7, xend = l*2/7, y = h*2 + u*2.5, yend = h*2 + u*2.5)) + 
    geom_segment(aes(x = l/7, xend = l*2/7, y = -h*2 + u*2.5, yend = -h*2 + u*2.5)) +
    geom_segment(aes(x = l*2/7, xend = l*3/7, y = h + u*2.5, yend = h + u*2.5)) + 
    geom_segment(aes(x = l*2/7, xend = l*3/7, y = -h + u*2.5, yend = -h + u*2.5)) +
    geom_segment(aes(x = l*2/7, xend = l*3/7, y = h*3 + u*2.5, yend = h*3 + u*2.5)) + 
    geom_segment(aes(x = l*2/7, xend = l*3/7, y = -h*3 + u*2.5, yend = -h*3 + u*2.5)) +
    geom_segment(aes(x = l/7, xend = l/7, y = -h*2 + u*2.5, yend = h*2 + u*2.5), linetype = "dashed", color = "gray") +
    geom_segment(aes(x = l*2/7, xend = l*2/7, y = -h*3 + u*2.5, yend = h*3 + u*2.5), linetype = "dashed", color = "gray") +
    ###
    theme(axis.title = element_blank(),
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(1, -1, 1, -1), "mm"),
          plot.title = element_text(hjust = 0.5, size = 8, margin = margin(b = -5))) +
    ggtitle("Config")
  
  g1 = ggplot(data, aes(x = config, y = coverage, fill = theta.type)) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) + 
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          plot.margin = unit(c(1, 0, 1, 0), "mm"),
          legend.position = "none", 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = .5, size = 8, margin = margin(b = -5)),
          panel.border = element_blank()) +
    geom_text(aes(x = config, y = coverage, label = coverage, group = theta.type), 
              hjust = 0, size = 2,
              position = position_dodge(width = 0.9)) +
    scale_y_reverse(limits = c(1, 0), breaks = c(0, 0.9)) + 
    scale_x_continuous(breaks = seq(0.5, 6.5, 1)) +
    scale_fill_manual(values = c(gray[3], gray[2], gray[1])) + 
    coord_flip() + 
    ggtitle("Coverage prob.")
  
  g2 <- ggplot(data = data, aes(x = config, y = ci.length, fill = theta.type)) + #xlab(NULL) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +# ggtitle("confidence interval length") +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1, 0, 1, -1), "mm"), 
          panel.grid.minor = element_blank(), 
          legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 8, margin = margin(b = -5)),
          panel.border = element_blank()) +
    geom_text(aes(x = config, y = ci.length, label = ci.length), 
              hjust = 1, size = 2,
              position = position_dodge(width = 0.9)) + 
    scale_y_continuous(limits = c(0, max(data$ci.length) * 1.01), breaks = c(0)) + 
    scale_fill_manual(values = c(gray[3], gray[2], gray[1])) + 
    scale_x_continuous(breaks = seq(0.5, 6.5, 1)) +
    coord_flip() + 
    ggtitle("CI length")
  
  gg1 <- ggplot_gtable(ggplot_build(g1))
  gg2 <- ggplot_gtable(ggplot_build(g2))
  gg.mid <- ggplot_gtable(ggplot_build(g.mid))
  
  grid.arrange(gg1, gg.mid, gg2, ncol = 3, widths = c(3/8, 2/8, 3/8))
  
  if (tikz) {
    tikz(file = filename, width = 3.1, height = 3)
    grid.arrange(gg1, gg.mid, gg2, ncol = 3, widths = c(3/8, 2/8, 3/8))
    dev.off()
  } else {
    grid.arrange(gg1, gg.mid, gg2, ncol = 3, widths = c(3/8, 2/8, 3/8))
  }
}

MakePlot(50, "lm", TRUE)
MakePlot(50, "logistic", TRUE)

# legend
tikz("./figs/legend_coverage.tex", width = 3, height = 0.15)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 1.1), ylim = c(0, 0.1))
polygon(c(0, 0.15, 0.15, 0), c(0, 0, 0.1, 0.1), col = gray[3], border = NA)
text(0.15, 0.05, label = 'null', pos = 4, cex = 0.7)
polygon(c(0.4, 0.55, 0.55, 0.4), c(0, 0, 0.1, 0.1), col = gray[2], border = NA)
text(0.55, 0.05, label = 'sparse', pos = 4, cex = 0.7)
polygon(c(0.8, 0.95, 0.95, 0.8), c(0, 0, 0.1, 0.1), col = gray[1], border = NA)
text(0.95, 0.05, label = 'dense', pos = 4, cex = 0.7)
dev.off()
