#
#  Function 'split_and_align' takes two input series, maps them
#  to the frequency domain, determines regions of high coherence,
#  then phase aligns the regions of x2 to closely match those of 
#  x1 in phase. The second series, x2, is then inverse-mapped back
#  to the time domain. The result should be an increased correlation
#  between the two series for the frequency region above the 
#  cutoff freq_cut.
#
"split_and_align" <- function(x1, x2, freq_cut, msc_sig, NW, K, dT) {
    N <- length(x1)
    stopifnot(length(x1) == length(x2), is.numeric(NW), is.numeric(K), is.numeric(dT))
    stopifnot(K <= 2 * NW, NW / N < 0.5)
     
    nFFT <- 2^(ceiling(log(N, 2)) + 1)
    x1 <- centre(x1, nw = NW, k = K, deltaT = dT)
    x2 <- centre(x2, nw = NW, k = K, deltaT = dT)
    s1 <- spec.mtm(x1, deltat = dT, Ftest = TRUE, returnInternals = TRUE, nw = NW, k = K,
                   plot = FALSE, nFFT = nFFT)
    s2 <- spec.mtm(x2, deltat = dT, Ftest = TRUE, returnInternals = TRUE, nw = NW, k = K,
                   plot = FALSE, nFFT = nFFT)
    coh <- mtm.coh(s1, s2, plot = FALSE)
    Wbin <- min(which(s1$freq >= 3 / s1$origin.n / s1$mtm$deltaT)) 
 
    #
    #  The filter for "which areas are important" is based on the coherence between
    #  the two series; ignoring Ftests and structure, since we're using complex-demod
    #  to grab the two "areas", approximately phase-align them, then put them back
    #
    peaks <- findPeaks(coh$NTmsc)  # binary "peak" or "not peak" label
    peaks2 <- filterPeaksW(peaks, coh$msc, Wbin)  # filter to remove twinned peaks within 2W bands
    while(sum(peaks2) < sum(peaks)) {
        peaks <- peaks2
        peaks2 <- filterPeaksW(peaks, coh$msc, Wbin)
    }
    peaks2 <- which(peaks2 == 1)
    mask1 <- peaks2[which(coh$msc[peaks2] > findMSCsig(coh, msc_sig)$msc)]

    #
    #  Now have length(mask1) frequencies which are central to (-W,W) bins which 
    #  are msc_sig (%) coherent or better. This is fairly vague, but hopefully good enough.
    #  These bands now are pulled out, one at a time (via complex demod) and compared
    #  to similarly pulled out response bands to determine the phase differences.
    #
    M <- length(coh$freq)
    egnC <- s2$mtm$eigenCoefs
    WbinUse <- floor(Wbin / 2)

    #  cut to 1 / freq_cut seconds period, only consider bands above this
    if(length(mask1) >= 1) {
      if(min(coh$freq[mask1]) < freq_cut) {
        jMin <- max(which(coh$freq[mask1] < freq_cut ))
      } else {
        jMin <- 0
      }
    } else {
      jMin <- 0
    }

    x2_zp <- c(x2, rep(0, nFFT - N))
    x2_pgrm <- fft(x2_zp)

    if(length(mask1) >= jMin + 1) {
      for(j in (jMin + 1):length(mask1)) {
        lowerB <- max(1, mask1[j] - WbinUse)
        upperB <- min(M, mask1[j] + WbinUse)
        rng <- lowerB:upperB
        phases <- coh$ph[rng] %% 360
        phases[phases > 179.9] <- phases[phases > 179.9] - 360
        phases <- phases / 180 * pi  # need phases in radians for correction
        offsetPh <- complex(real = cos(phases), imaginary = sin(phases))

        # start with x2, flip it into freq domain with regular periodogram,
        # then apply this phase. Rinse-repeat over each of the adjustments.
        x2_pgrm[rng] <- x2_pgrm[rng] * offsetPh
      }

      # Now have phase-adjusted array, continue and invert to "new" series
      x2_pgrm[(nFFT/2+2):nFFT] <- Conj(x2_pgrm[(nFFT/2):2])
      flip_back <- Re(fft(x2_pgrm, inverse = TRUE))[1:N] / nFFT
    } else {
      # didn't find any coherent peaks to phase-align; just return x2
      return(x2)
    }
    return(flip_back)
}



