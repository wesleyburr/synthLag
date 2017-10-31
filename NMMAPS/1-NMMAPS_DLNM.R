#
#  Using the AHI database, compute all-year (1987-2000, to match NMMAPS)
#  risks, lags 0, 1 and 2. Then compute the synthetic lag versions of this.
#  Then compute DLNMs.
#

library("AHItools")
library("dlnm")
library("splines")

max_iter <- 100
epsilon <- 1e-9
dfs <- 12
K <- 14 * 12
N <- 5114
W <- round(K / N * 365.2425) / 2  # round to nearest 0.5
timeBasis_slp <- load_time_basis(W = W, K = K, N = N)
gam.stmt <- c("glm(as.formula(",
              paste0("), data = dat, family = quasi(log, mu),
              na.action = na.exclude, epsilon = ", epsilon,
              ", maxit = ", max_iter, ")"))

################################################################################
#
#  NMMAPS using DLNM (100 biggest cities, not all available)
#
################################################################################

stop("You must have the NMMAPS city-specific complete databases stored in the NMMAPS folder. They are 3.4GB in size, and 
    are not distributed with this code.")
all_cities <- list.files("./NMMAPS/")
all_cities <- all_cities[-which(all_cities == "cities.rda")]

obj_names <- unlist(lapply(strsplit(all_cities, "\\."), "[[", 1))

res_NMMAPS <- data.frame(allRRlow  = rep(NA, length(obj_names)),
                         allRRfit  = rep(NA, length(obj_names)),
                         allRRhigh = rep(NA, length(obj_names)))
load("./NMMAPS/cities.rda")
row.names(res_NMMAPS) <- cities$city

model1.stmt <- "cvd ~ model.o3 + tmpd + timeBasis_slp + as.factor(dow)"

for(j in 1:length(all_cities)) {
  load(paste0("./NMMAPS/", all_cities[j]))
  dat <- get(obj_names[j])
  dat2 <- cbind(as.numeric(substr(dat$date[1:N], 1, 4)), dat[1:N, c("cvd", "tmpd", "dow", "o3mean")])
  names(dat2)[1] <- "year"
  dat2$cvd <- dat[1:N, "cvd"] + dat[(N+1):(2*N), "cvd"] + dat[(2*N+1):(3*N), "cvd"]
  dat <- dat2
  if(length(which(!is.na(dat$o3mean))) > 0) {
    model.o3 <- crossbasis(dat$o3mean, lag = 15,
                           argvar = list(fun = "lin"),
                           arglag = list(fun = "poly", degree = 4))

    fit <- eval(parse(text = paste0(gam.stmt[1], "model1.stmt", gam.stmt[2])))
    pred1.o3 <- crosspred(model.o3, fit, at = 0:20, bylag = 0.2,
                          cumul = TRUE)

    res_NMMAPS[j, ] <- c(pred1.o3$allRRlow[2], pred1.o3$allRRfit[2], pred1.o3$allRRhigh[2])
  }
}
names(res_NMMAPS) <- names(res_AHI) <- c("allRRlow", "allRRfit", "allRRhigh")

################################################################################
#
#  Next, NMMAPS using synth-lag - manually
#
################################################################################

res_NMMAPS_synth <- data.frame(allRRlow  = rep(NA, length(obj_names)),
                               allRRfit  = rep(NA, length(obj_names)),
                               allRRhigh = rep(NA, length(obj_names)))
load("./NMMAPS/cities.rda")
row.names(res_NMMAPS) <- cities$city

model1.stmt <- "cvd ~ o3mean_synth + tmpd + timeBasis_slp + as.factor(dow)"
synth_sig <- 0.95
synth_cut = 1/60/86400
res_NMMAPS_synth <- data.frame(allRRlow  = rep(NA, length(obj_names)),
                               allRRfit  = rep(NA, length(obj_names)),
                               allRRhigh = rep(NA, length(obj_names)))

for(j in 1:length(all_cities)) {
  load(paste0("./NMMAPS/", all_cities[j]))
  dat <- get(obj_names[j])
  dat2 <- cbind(as.numeric(substr(dat$date[1:N], 1, 4)), dat[1:N, c("cvd", "tmpd", "dow", "o3mean")])
  names(dat2)[1] <- "year"
  dat2$cvd <- dat[1:N, "cvd"] + dat[(N+1):(2*N), "cvd"] + dat[(2*N+1):(3*N), "cvd"]
  dat <- dat2

  if(length(which(is.na(dat$o3mean))) == 0) {
    dat$cvd[dat$cvd == 0] <- 0.01
    new_synth_lag_series <- matrix(data = NA, nrow = dim(dat)[1],
                                          ncol = 1)

    # annual estimates: use the synth_span parameter
    synth_span <- 5
    for(yr in 1987:2000) {
     data_tmp <- dat[dat$year >= yr - synth_span / 2 & dat$year <= yr + synth_span / 2, ]
     new_synth_lag_series[dat$year == yr,  1] <-
              split_and_align(log(data_tmp$cvd), data_tmp$o3mean, freq_cut = synth_cut,
                              msc_sig = synth_sig, NW = 3, K = 5, dT = 86400)[data_tmp$year == yr]
    }
    dat <- cbind(dat, new_synth_lag_series)
    names(dat)[length(names(dat))] <- "o3mean_synth"

    fit <- eval(parse(text = paste0(gam.stmt[1], "model1.stmt", gam.stmt[2])))

    coefs <- summary.glm(fit)$coefficients["o3mean_synth", ]
    z <- qnorm(1-(1-0.95)/2)
    res_NMMAPS_synth[j, ] <- c(coefs[1] - coefs[2] * z, coefs[1], coefs[1] + coefs[2] * z)
  }
}
names(res_NMMAPS_synth) <- c("allRRlow", "allRRfit", "allRRhigh")
row.names(res_NMMAPS_synth) <- cities$city

res_NMMAPS <- res_NMMAPS[!is.na(res_NMMAPS_synth[, 1]), ]
res_NMMAPS_synth <- res_NMMAPS_synth[!is.na(res_NMMAPS_synth[, 1]), ]

################################################################################
#
#  Results done: set of 22 AHI dlnm & synth, set of NMMAPS (21) dlnm & synth
#
################################################################################
save(res_NMMAPS, file = "res_NMMAPS_dlnm.RDa") 
save(res_NMMAPS_synth, file = "res_NMMAPS_synth.RDa")

