#
#  Simulation 3
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
phi_var1 <- seq(-pi, pi, pi/10)
cors <- vector(length = M)
lmod <- vector(length = M)

# Estimate linear model, no modification
comp <- function(X, M) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 0.5)
    wn2 <- rnorm(N, 0, 0.5)
    x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + X) 
    x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) 
    lmod[j] <- lm(x2 ~ x1)$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}

# Estimate linear model, lagging by-day
comp_lag1 <- function(X, M) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 0.5)
    wn2 <- rnorm(N, 0, 0.5)
    x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + X) 
    x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) 
    lmod[j] <- lm(x2[-N] ~ x1[-1])$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}

comp_lag2 <- function(X, M) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 0.5)
    wn2 <- rnorm(N, 0, 0.5)
    x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + X) 
    x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) 
    lmod[j] <- lm(x2[-c(N, N-1)] ~ x1[-c(1, 2)])$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}

comp_lag3 <- function(X, M) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 0.5)
    wn2 <- rnorm(N, 0, 0.5)
    x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + X) 
    x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) 
    lmod[j] <- lm(x2[-c(N:(N-2))] ~ x1[-c(1:3)])$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}

# Estimate linear model, synthetic lag
comp_syn <- function(X, M, synth_cut = 1e-10, synth_sig = 0.999) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 0.5)
    wn2 <- rnorm(N, 0, 0.5)
    x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + X) 
    x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) 
    x1s <- split_and_align(x2, x1, freq_cut = synth_cut,
                           msc_sig = synth_sig, NW = 3, K = 5, 
                           dT = 86400
                          )
    lmod[j] <- lm(x2 ~ x1s)$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}

#
#  This is not fast, it's brute-force. There's many ways to improve the speed, but it's
#  a demonstration.
#
res <- rep(NA, length(phi_var1))

results <- data.frame(phi = phi_var1, 
                      coef0 = res, coef1 = res, coef2 = res, coef3 = res, coefs = res,
                      ymin0 = res, ymin1 = res, ymin2 = res, ymin3 = res, ymins = res,
                      ymax0 = res, ymax1 = res, ymax2 = res, ymax3 = res, ymaxs = res)
    
for(j in 1:length(phi_var1)) {
  vals <- comp(phi_var1[j], M)
  vals_syn <- comp_syn(phi_var1[j], M)
  vals_lag1 <- comp_lag1(phi_var1[j], M)
  vals_lag2 <- comp_lag2(phi_var1[j], M)
  vals_lag3 <- comp_lag3(phi_var1[j], M)

  results[j, ] <- c(phi_var1[j], vals[1], vals_lag1[1], vals_lag2[1], vals_lag3[1], vals_syn[1],
                    vals[1] + qnorm(0.025) * vals[2],
                    vals_lag1[1] + qnorm(0.025) * vals_lag1[2],
                    vals_lag2[1] + qnorm(0.025) * vals_lag2[2],
                    vals_lag3[1] + qnorm(0.025) * vals_lag3[2],
                    vals_syn[1] + qnorm(0.025) * vals_syn[2],
                    vals[1] + qnorm(0.975) * vals[2],
                    vals_lag1[1] + qnorm(0.975) * vals_lag1[2],
                    vals_lag2[1] + qnorm(0.975) * vals_lag2[2],
                    vals_lag3[1] + qnorm(0.975) * vals_lag3[2],
                    vals_syn[1] + qnorm(0.975) * vals_syn[2])
}

p <- ggplot(results, aes(phi)) + 
       geom_ribbon(aes(ymin = ymins, ymax = ymaxs), fill = "grey80") + 
       geom_hline(yintercept = 0, lty = 1) + 
       geom_line(aes(y = coef0), lty = 1) + 
       geom_line(aes(y = coef1), lty = 2) + 
       geom_line(aes(y = coef2), lty = 3) + 
       geom_line(aes(y = coef3), lty = 4) + 
       geom_line(aes(y = coefs), lwd = 2) + 
       ylab("Coefficient of Linear Model") +
       xlab(expression(phi))
plot(p)

pdf(file = "sim3.pdf", width = 8, height = 6)
par(mar = c(4, 4, 0.5, 0.5))
plot(p)
dev.off()


