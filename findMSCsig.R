#
#  findMSCsig : function to find the MSC which corresponds to a significance
#

"findMSCsig" <- function(cohObj, percentG) {
  stopifnot(class(cohObj) == "mtm.coh")

  k <- cohObj$mtm$k
  trnrm_ <- sqrt(2 * k - 2)
  MSqCoh <- 1.0 - (1.0 - percentG)**(1.0 / as.double(k - 1))

  NTmsc <- trnrm_ * log( (1.0 + sqrt(MSqCoh)) / (1.0 - sqrt(MSqCoh)) ) / 2.0
  msc <- MSqCoh
  obj <- list(NTmsc, msc, percentG)
  names(obj) <- c("NTmsc", "msc", "percentG")
  return(obj)
}


