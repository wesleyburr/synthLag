#
#  findPeaks : locate local min/max in a time series
#
"findPeaks" <- function(dat) {
  N <- length(dat)
  stopifnot(is.numeric(dat))

  locMax <- rep(0, N)

  for(j in 2:(N-1)) {
    if(dat[j] > dat[j-1] & dat[j] > dat[j+1]) {
      locMax[j] <- 1
    }
  }
  locMax
}


