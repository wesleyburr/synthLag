# 
#  explore.R
#  * issues with variability of phase-aligned estimates at -90 to -180 and +90 to +180 offsets
#
library("ggplot2")
library("reshape2")
library("multitaper")
source("findPeaks.R")
source("filterPeaksW.R")
source("findMSCsig.R")

N <- 1000
tarr <- 1:N
phi <- 5/4 * pi / 2

wn1 <- rnorm(N, 0, 1)
wn2 <- rnorm(N, 0, 1)
x1 <- wn1 + 2 * cos(2 * pi * 1/4 * tarr + 0.0) 
x2 <- wn2 + 2 * cos(2 * pi * 1/4 * tarr + phi) 

# arguments to split_and_align
NW <- 4
K <- 7
dT <- 86400
msc_sig <- 0.999

N <- length(x1)
x1 <- centre(x1, nw = NW, k = K, deltaT = dT)
x2 <- centre(x2, nw = NW, k = K, deltaT = dT)
s1 <- spec.mtm(x1, deltat = dT, Ftest = TRUE, returnInternals = TRUE, nw = NW, k = K,
               plot = FALSE, adaptiveWeighting = FALSE)
s2 <- spec.mtm(x2, deltat = dT, Ftest = TRUE, returnInternals = TRUE, nw = NW, k = K,
               plot = FALSE, adaptiveWeighting = FALSE)
coh <- mtm.coh(s1, s2, plot = FALSE)
Wbin <- min(which(s1$freq >= 3 / s1$origin.n / s1$mtm$deltaT)) 

peaks <- findPeaks(coh$NTmsc)  # binary "peak" or "not peak" label
peaks2 <- filterPeaksW(peaks, coh$msc, Wbin)  # filter to remove twinned peaks within 2W bands
while(sum(peaks2) < sum(peaks)) {
    peaks <- peaks2
    peaks2 <- filterPeaksW(peaks, coh$msc, Wbin)
}
peaks2 <- which(peaks2 == 1)
mask1 <- peaks2[which(coh$msc[peaks2] > findMSCsig(coh, msc_sig)$msc)]

M <- length(coh$freq)
egnC <- s2$mtm$eigenCoefs
WbinUse <- floor(Wbin / 2)

for(j in 1:length(mask1)) {
  lowerB <- max(1, mask1[j] - WbinUse)
  upperB <- min(M, mask1[j] + WbinUse)
  rng <- lowerB:upperB
  phases <- coh$ph[rng] %% 360
  phases[phases > 179.9] <- phases[phases > 179.9] - 360
  phases <- phases / 180 * pi  # need phases in radians for correction
  offsetPh <- complex(real = cos(phases), imaginary = sin(phases))
  egnC[rng, ] <- apply(egnC[rng, ], MAR = 2, FUN = function(x) { x * offsetPh })
}

# make full CFT array
cft <- rbind(egnC, Conj(egnC[(s2$mtm$nfreqs-1):2, ]))
# invert
inv <- apply(cft, MAR = 2, FUN = function(x) { 1 / s2$mtm$nFFT * Re(fft(x, inverse = TRUE))[1:N] })

#  Invert by using the trick of the Ftest: \sum_{k} x * v_k^2 / \sum_{k} v_k^2 
scaleFactor <- s2$mtm$dpss$v * sqrt(s2$mtm$deltaT)
fixedX2 <- rowSums(inv * scaleFactor) / rowSums(scaleFactor * scaleFactor)

#
#  Why the fuck is the amplitude of the peak decreasing like this? It makes no sense!
#

#
#  Try the multiplication on the raw form: just using FFT, and rotate
#
nFFT <- 2048
test <- x2
test_zp <- c(x2, rep(0, nFFT - 1000))
sp_test <- fft(test_zp)
# try to do the rotation
sp_test[rng] <- sp_test[rng] * offsetPh
sp_test[(nFFT/2+2):nFFT] <- Conj(sp_test[(nFFT/2):2])
test_back <- Re(fft(sp_test, inverse = TRUE))[1:1000] / nFFT
plot(test - test_back, type = "l")

#
#  This isn't perfect: you lose a little bit of the magnitude, but not as much as
#  happens when you try to do it through multitaper
#
