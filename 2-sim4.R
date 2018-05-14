#
#  Simulation 4
#  * independent WN(0, \sigma^2) base
#  * 10 periodic component, amplitudes varying from 0.5 to 3 (randomly)
#  * equal \omega terms across the 10 periodic components
#  * \phi terms vary randomly from [-pi, pi]
#  * true 'relationship' estimated by using 0-phase offsets
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
M <- 250
tarr <- 1:N

freqs <- 1 / sample(1:100, 10, replace = FALSE)
amps <- sample(seq(0.5, 3.0, 0.1), 10)
synth_cut <- 1e-10
synth_sig = 0.90

lmods0 <- lmods1 <- lmods2 <- lmods3 <- lmods_syn <- vector(length = M)

# iterate through realizations
for(j in 1:M) {
  x1 <- x1_noPh <- rnorm(N, 0, 1)
  x2 <- rnorm(N, 0, 1)
  for(m in 1:length(freqs)) {
    x1 <- x1 + amps[m] * cos(2 * pi * freqs[m] * tarr + runif(1, min = -pi, max = pi))
    x1_noPh <- x1_noPh + amps[m] * cos(2 * pi * freqs[m] * tarr + 0.0)
    x2 <- x2 + amps[m] * cos(2 * pi * freqs[m] * tarr + 0.0)
  }
  true_lmod <- lm(x2 ~ x1_noPh)$coefficients[2]

  # Integer-valued lags
  lmods0[j] <- lm(x2 ~ x1)$coefficients[2] / true_lmod * 100   # record the percent of true, +/-
  lmods1[j] <- lm(x2[-N] ~ x1[-1])$coefficients[2] / true_lmod * 100   # record the percent of true, +/-
  lmods2[j] <- lm(x2[-c(N, N-1)] ~ x1[-c(1, 2)])$coefficients[2] / true_lmod * 100   # record the percent of true, +/-
  lmods3[j] <- lm(x2[-(N:(N-2))] ~ x1[-(1:3)])$coefficients[2] / true_lmod * 100   # record the percent of true, +/-

  # Synthetic lag
  x1s <- split_and_align(x2, x1, freq_cut = synth_cut,
                         msc_sig = synth_sig, NW = 3, K = 5, 
                         dT = 86400
                        )
  lmods_syn[j] <- lm(x2 ~ x1s)$coefficients[2] / true_lmod * 100
}

res <- matrix(c(quantile(lmods0, probs = c(0.025, 0.50, 0.975)), 
                quantile(lmods1, probs = c(0.025, 0.50, 0.975)), 
                quantile(lmods2, probs = c(0.025, 0.50, 0.975)), 
                quantile(lmods3, probs = c(0.025, 0.50, 0.975)), 
                quantile(lmods_syn, probs = c(0.025, 0.50, 0.975))),
              byrow = TRUE, nrow = 5, ncol = 3)
res <- as.data.frame(res)
rownames(res) <- c("Lag-0", "Lag-1", "Lag-2", "Lag-3", "SynthLag")
colnames(res) <- c("2.5%", "50%", "97.5%")

# We actually have the simulation-by-simulation results as percentages, so take
# the synthLag result by case, and subtract the other 3 from it, and report this as
# additional columns

res2 <- res
res2[1, ] <- quantile(lmods_syn - lmods0, probs = c(0.025, 0.50, 0.975))
res2[2, ] <- quantile(lmods_syn - lmods1, probs = c(0.025, 0.50, 0.975))
res2[3, ] <- quantile(lmods_syn - lmods2, probs = c(0.025, 0.50, 0.975))
res2[4, ] <- quantile(lmods_syn - lmods3, probs = c(0.025, 0.50, 0.975))
res2[5, ] <- c(NA, NA, NA)



