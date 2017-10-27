#
#  Filter a peaks list to find local maximums
#
"filterPeaksW" <- function(peaks, orig, W) {

  totJ <- 0
  stopifnot(is.numeric(peaks), max(peaks) == 1, min(peaks) == 0,
            length(peaks) == length(orig))
  N <- length(peaks)
  peaks2 <- peaks

  j <- W+1
  while(j < (N-W)) {
      if(sum( peaks[(j-W):(j+W)] ) > 1) {
        peaks2[(j-W):(j+W)] <- rep(0, 2 * W + 1)      
        actualMax <- which(orig[(j-W):(j+W)] == max(orig[(j-W):(j+W)]))
        peaks2[j-W - 1 + actualMax] <- 1
        j <- j + (2 * W + 1)
      } else {
          j <- j + 1
      }
  }
  peaks2
}


