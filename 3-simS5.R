#
#  Simulation S4
#  * independent WN(0, \sigma^2) base
#  * NO periodic components
#  * the "null" case
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
lmod_s <- lmod

#
#  This is not fast, it's brute-force. There's many ways to improve the speed, but it's
#  a demonstration.
#

# Estimate linear model, synthetic lag + default
M <- 200
synth_cut = 1e-10
synth_sig = 0.999

for(j in 1:M) {
  wn1 <- rnorm(N, 0, 0.5)
  wn2 <- rnorm(N, 0, 0.5)
  x1 <- wn1
  x2 <- wn2
  x1s <- split_and_align(x2, x1, freq_cut = synth_cut,
                         msc_sig = synth_sig, NW = 3, K = 5, 
                         dT = 86400
                        )
  lmod[j] <- lm(x2 ~ x1)$coefficients[2]
  lmod_s[j] <- lm(x2 ~ x1s)$coefficients[2]
}
print(sort(lmod)[100])
print(sort(lmod)[5])
print(sort(lmod)[195])

print(sort(lmod_s)[100])
print(sort(lmod_s)[5])
print(sort(lmod_s)[195])


N <- 5114
M <- 200
tarr <- 1:N

for(j in 1:M) {
  wn1 <- rnorm(N, 0, 0.5)
  wn2 <- rnorm(N, 0, 0.5)
  x1 <- wn1
  x2 <- wn2
  x1s <- split_and_align(x2, x1, freq_cut = synth_cut,
                         msc_sig = synth_sig, NW = 3, K = 5, 
                         dT = 86400
                        )
  lmod[j] <- lm(x2 ~ x1)$coefficients[2]
  lmod_s[j] <- lm(x2 ~ x1s)$coefficients[2]
}

print(sort(lmod)[100])
print(sort(lmod)[5])
print(sort(lmod)[195])

print(sort(lmod_s)[100])
print(sort(lmod_s)[5])
print(sort(lmod_s)[195])


