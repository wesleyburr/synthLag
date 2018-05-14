#
#  Simulation 4
#  * independent WN(0, \sigma^2) base
#  * two periodic component
#  * equal \omega terms
#  * \phi_1 and \phi_2 for predictor vary in discrete steps seq(-pi, pi, pi/6)
#
library("ggplot2")
library("reshape2")
library("multitaper")
source("split_and_align.R")
source("findPeaks.R")
source("filterPeaksW.R")
source("findMSCsig.R")

set.seed(67)

N <- 1000
M <- 100
tarr <- 1:N
phi_var1 <- seq(-pi, pi, pi/6)
phi_var2 <- seq(-pi, pi, pi/6)
cors <- vector(length = M)
lmod <- vector(length = M)

# Estimate linear model, no modification
comp <- function(X, Y, M) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 1)
    wn2 <- rnorm(N, 0, 1)
    x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + X) + 2 * cos(2 * pi * 1/6 * tarr + Y)
    x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) + 2 * cos(2 * pi * 1/6 * tarr + 0.0)
    lmod[j] <- lm(x2 ~ x1)$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}

# Estimate linear model, synthetic lag
comp_syn <- function(X, Y, M, synth_cut = 1e-10, synth_sig = 0.999) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 1)
    wn2 <- rnorm(N, 0, 1)
    x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + X) + 2 * cos(2 * pi * 1/6 * tarr + Y)
    x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) + 2 * cos(2 * pi * 1/6 * tarr + 0.0)
    x1s <- split_and_align(x2, x1, freq_cut = synth_cut,
                           msc_sig = synth_sig, NW = 3, K = 5, 
                           dT = 86400
                          )
    lmod[j] <- lm(x2 ~ x1s)$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}

res1 <- res2 <- res1s <- res2s <- matrix(data = NA, nrow = length(phi_var1), ncol = length(phi_var2))
#
#  This is not fast, it's brute-force. There's many ways to improve the speed, but it's
#  a demonstration.
#
for(j in 1:length(phi_var1)) {
  cat(".")
  for(k in 1:length(phi_var2)) {
    vals <- comp(phi_var1[j], phi_var2[k], M)
    vals_syn <- comp_syn(phi_var1[j], phi_var2[k], M)

    res1[j, k] <- vals[1]
    res2[j, k] <- vals[2]
    res1s[j, k] <- vals_syn[1]
    res2s[j, k] <- vals_syn[2]
  }
}
cat("\n")

colnames(res1) <- colnames(res1s) <- phi_var1
rownames(res1) <- rownames(res1s) <- phi_var2
melted_lmod1 <- melt(res1)
melted_lmod1s <- melt(res1s)

melted_lmod1s_perc <- melted_lmod1s
melted_lmod1s_perc$value <- melted_lmod1s_perc$value / 0.80

q <- ggplot(data = melted_lmod1s_perc, aes(x=Var1, y=Var2, fill=value)) + 
     geom_tile(aes(fill = value), colour = "white") +
     scale_fill_gradient(low = "white", high = "steelblue") +
     ylab(expression(phi[2])) +
     xlab(expression(phi[1]))


