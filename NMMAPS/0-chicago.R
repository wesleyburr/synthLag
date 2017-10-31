multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#
#  Analyze Chicago NMMAPS data using synthetic lag
#
library("lubridate")
library("ggplot2")
load("./NMMAPS/chic.rda")
chic <- chic[, c("date", "dow", "agecat", "cvd", "death", "resp", "tmpd", "tmean",
                 "pm10tmean", "o3tmean", 
                 "l1pm10tmean", "l1o3tmean",
                 "l2pm10tmean", "l2o3tmean",
                 "l3pm10tmean", "l3o3tmean")]
# Split the date
Time  <- rep(1:(dim(chic)[1] / 3), 3)
Year  <- as.numeric(substr(chic$date, 1, 4))
Month <- as.numeric(substr(chic$date, 5, 6))
Day   <- as.numeric(substr(chic$date, 7, 8))
chic <- data.frame(chic, Time, Year, Month, Day)

# Compress age categories
n_rec <- dim(chic)[1] / 3
data.spec <- data.frame(chic[1:n_rec, c("date", "dow", "tmpd", "tmean", "pm10tmean", "o3tmean", "l1pm10tmean", "l1o3tmean",
                                        "l2pm10tmean", "l2o3tmean", "l3pm10tmean", "l3o3tmean", "Time", "Year", "Month", "Day")],
                        cvd = chic[1:n_rec, "cvd"] + chic[(n_rec+1):(2*n_rec), "cvd"] + chic[(n_rec*2+1):(3*n_rec), "cvd"],
                        death = chic[1:n_rec, "death"] + chic[(n_rec+1):(2*n_rec), "death"] + chic[(n_rec*2+1):(3*n_rec), "death"],
                        resp = chic[1:n_rec, "resp"] + chic[(n_rec+1):(2*n_rec), "resp"] + chic[(n_rec*2+1):(3*n_rec), "resp"])

#
#  Function to compute simple model coefficients
#
simple_compute <- function(data.spec, dfs, epsilon, bf.epsilon, 
                           max_iter = 100, bf.max_iter = 50,
                           resp, airpoll, tempVar) {

    require("gam")
    require("AHItools")
    N <- dim(data.spec)[1]
    timeBasis_ns <- ns(data.spec[, "Time"], df = dfs * length(unique(data.spec[, "Year"]))) 
    col_vect <- which(apply(timeBasis_ns, MAR = 2, FUN = max) > 1e-3)
    timeBasis_ns <- timeBasis_ns[, col_vect]

    # Issue: to match up the time bases, we want K to be the same ... but for a good SLP basis, we should
    # be using 2*N*W (approx) basis vectors: W = K / 2 / N
    K <- length(unique(data.spec[, "Year"])) * dfs
    W <- round(K / N * 365.2425) / 2  # round to nearest 0.5
    timeBasis_slp <- load_time_basis(W = W, K = K, N = N)
    col_vect <- which(apply(timeBasis_slp, MAR = 2, FUN = max) > 1e-3)
    timeBasis_slp <- timeBasis_slp[, col_vect]

    gam.stmt <- c("gam(as.formula(", paste0("), data = data.spec, family = quasi(log, mu), 
                      na.action = na.exclude, epsilon = ", epsilon, ", bf.epsilon = ", bf.epsilon, ", maxit = ", max_iter,
                      ", bf.maxit = ", bf.max_iter, ")"))
    model.stmt1 <- paste0(resp, " ~ ", airpoll, " + as.factor(dow) + ",
                          tempVar, " + timeBasis_ns")
    model.stmt2 <- paste0(resp, " ~ ", airpoll, " + as.factor(dow) + ",
                          tempVar, " + timeBasis_slp")

    fit1 <- tryCatch({
       eval(parse(text = paste0(gam.stmt[1], "model.stmt1", gam.stmt[2])))
    }, warning = function(w) {
      NULL 
    }, error = function(e) {
      NULL 
    }, finally = {
      NULL 
    })

    fit2 <- tryCatch({
       eval(parse(text = paste0(gam.stmt[1], "model.stmt2", gam.stmt[2])))
    }, warning = function(w) {
      NULL 
    }, error = function(e) {
      NULL 
    }, finally = {
      NULL 
    })
    two <- rep(NA, 2)
    returned <- data.frame(Estimate = two, 
                           Std.Dev = two,
                           pval = two,
                           iter = two)
    returned[1, "Estimate"] <- summary.glm(fit1)$coefficients[airpoll, "Estimate"]
    returned[1, "Std.Dev"] <- summary.glm(fit1)$coefficients[airpoll, "Std. Error"]
    returned[1, "pval"] <- summary.glm(fit1)$coefficients[airpoll, "Pr(>|t|)"]
    returned[1, "iter"] <- fit1$iter 

    returned[2, "Estimate"] <- summary.glm(fit2)$coefficients[airpoll, "Estimate"]
    returned[2, "Std.Dev"] <- summary.glm(fit2)$coefficients[airpoll, "Std. Error"]
    returned[2, "pval"] <- summary.glm(fit2)$coefficients[airpoll, "Pr(>|t|)"]
    returned[2, "iter"] <- fit2$iter 
    
    returned
 }

results_classic <- vector("list", 8)
results_classic[[1]] <- simple_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                       max_iter = 100, bf.max_iter = 100, 
                                       resp = "death", airpoll = "pm10tmean", tempVar = "tmpd")
results_classic[[2]] <- simple_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                       max_iter = 100, bf.max_iter = 100, 
                                       resp = "death", airpoll = "l1pm10tmean", tempVar = "tmpd")
results_classic[[3]] <- simple_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                       max_iter = 100, bf.max_iter = 100, 
                                       resp = "death", airpoll = "l2pm10tmean", tempVar = "tmpd")
results_classic[[4]] <- simple_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                       max_iter = 100, bf.max_iter = 100, 
                                       resp = "death", airpoll = "l3pm10tmean", tempVar = "tmpd")

results_classic[[5]] <- simple_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                     max_iter = 100, bf.max_iter = 100, 
                                     resp = "death", airpoll = "pm10tmean", tempVar = "tmpd")
results_classic[[6]] <- simple_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                       max_iter = 100, bf.max_iter = 100, 
                                       resp = "death", airpoll = "l1pm10tmean", tempVar = "tmpd")
results_classic[[7]] <- simple_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                       max_iter = 100, bf.max_iter = 100, 
                                       resp = "death", airpoll = "l2pm10tmean", tempVar = "tmpd")
results_classic[[8]] <- simple_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                       max_iter = 100, bf.max_iter = 100, 
                                       resp = "death", airpoll = "l3pm10tmean", tempVar = "tmpd")


synth_compute <- function(data.spec, dfs, epsilon, bf.epsilon, 
                          max_iter = 100, bf.max_iter = 50,
                          resp, airpoll, tempVar, synth_cut = 1 / 60 / 86400,
                          synth_sig = 0.90) {

    require("gam")
    require("AHItools")
 
    # Create the synthetically lagged series
    data.spec[which(is.na(data.spec[, airpoll])), airpoll] <- mean(data.spec[, airpoll], na.rm = TRUE)
    new_series <- split_and_align(log(data.spec[, resp]), data.spec[, airpoll], freq_cut = synth_cut,
                                  msc_sig = synth_sig, NW = 3, K = 5, dT = 86400)
 
    data.spec <- data.frame(data.spec, airpoll_synth = new_series)

    N <- dim(data.spec)[1]
    timeBasis_ns <- ns(data.spec[, "Time"], df = dfs * length(unique(data.spec[, "Year"]))) 
    col_vect <- which(apply(timeBasis_ns, MAR = 2, FUN = max) > 1e-3)
    timeBasis_ns <- timeBasis_ns[, col_vect]

    # Issue: to match up the time bases, we want K to be the same ... but for a good SLP basis, we should
    # be using 2*N*W (approx) basis vectors: W = K / 2 / N
    K <- length(unique(data.spec[, "Year"])) * dfs
    W <- round(K / N * 365.2425) / 2  # round to nearest 0.5
    timeBasis_slp <- load_time_basis(W = W, K = K, N = N)
    col_vect <- which(apply(timeBasis_slp, MAR = 2, FUN = max) > 1e-3)
    timeBasis_slp <- timeBasis_slp[, col_vect]

    gam.stmt <- c("gam(as.formula(", paste0("), data = data.spec, family = quasi(log, mu), 
                      na.action = na.exclude, epsilon = ", epsilon, ", bf.epsilon = ", bf.epsilon, ", maxit = ", max_iter,
                      ", bf.maxit = ", bf.max_iter, ")"))
    model.stmt1 <- paste0(resp, " ~ airpoll_synth + as.factor(dow) + ",
                          tempVar, " + timeBasis_ns")
    model.stmt2 <- paste0(resp, " ~ airpoll_synth + as.factor(dow) + ",
                          tempVar, " + timeBasis_slp")

    fit1 <- tryCatch({
       eval(parse(text = paste0(gam.stmt[1], "model.stmt1", gam.stmt[2])))
    }, warning = function(w) {
      NULL 
    }, error = function(e) {
      NULL 
    }, finally = {
      NULL 
    })

    fit2 <- tryCatch({
       eval(parse(text = paste0(gam.stmt[1], "model.stmt2", gam.stmt[2])))
    }, warning = function(w) {
      NULL 
    }, error = function(e) {
      NULL 
    }, finally = {
      NULL 
    })
    two <- rep(NA, 2)
    returned <- data.frame(Estimate = two, 
                           Std.Dev = two,
                           pval = two,
                           iter = two)
    returned[1, "Estimate"] <- summary.glm(fit1)$coefficients["airpoll_synth", "Estimate"]
    returned[1, "Std.Dev"] <- summary.glm(fit1)$coefficients["airpoll_synth", "Std. Error"]
    returned[1, "pval"] <- summary.glm(fit1)$coefficients["airpoll_synth", "Pr(>|t|)"]
    returned[1, "iter"] <- fit1$iter 

    returned[2, "Estimate"] <- summary.glm(fit2)$coefficients["airpoll_synth", "Estimate"]
    returned[2, "Std.Dev"] <- summary.glm(fit2)$coefficients["airpoll_synth", "Std. Error"]
    returned[2, "pval"] <- summary.glm(fit2)$coefficients["airpoll_synth", "Pr(>|t|)"]
    returned[2, "iter"] <- fit2$iter 
    
    returned
}

results_synth <- vector("list", 8)
results_synth[[1]] <- synth_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                    max_iter = 100, bf.max_iter = 100, 
                                    resp = "death", airpoll = "pm10tmean", tempVar = "tmpd")
results_synth[[2]] <- synth_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                    max_iter = 100, bf.max_iter = 100, 
                                    resp = "death", airpoll = "l1pm10tmean", tempVar = "tmpd")
results_synth[[3]] <- synth_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                    max_iter = 100, bf.max_iter = 100, 
                                    resp = "death", airpoll = "l2pm10tmean", tempVar = "tmpd")
results_synth[[4]] <- synth_compute(data.spec = data.spec, dfs = 8, epsilon = 1e-7, bf.epsilon = 1e-7,
                                    max_iter = 100, bf.max_iter = 100, 
                                    resp = "death", airpoll = "l3pm10tmean", tempVar = "tmpd")

results_synth[[5]] <- synth_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                  max_iter = 100, bf.max_iter = 100, 
                                  resp = "death", airpoll = "pm10tmean", tempVar = "tmpd")
results_synth[[6]] <- synth_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                    max_iter = 100, bf.max_iter = 100, 
                                    resp = "death", airpoll = "l1pm10tmean", tempVar = "tmpd")
results_synth[[7]] <- synth_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                    max_iter = 100, bf.max_iter = 100, 
                                    resp = "death", airpoll = "l2pm10tmean", tempVar = "tmpd")
results_synth[[8]] <- synth_compute(data.spec = data.spec, dfs = 12, epsilon = 1e-7, bf.epsilon = 1e-7,
                                    max_iter = 100, bf.max_iter = 100, 
                                    resp = "death", airpoll = "l3pm10tmean", tempVar = "tmpd")

#
#  points of interest for the talk
#
pm10_ci <- matrix(data = NA, nrow = 4, ncol = 2)
pm10_ci[1, ] <- unlist(results_classic[[1]][1, 1:2])
pm10_ci[2, ] <- unlist(results_classic[[2]][1, 1:2])
pm10_ci[3, ] <- unlist(results_classic[[3]][1, 1:2])
pm10_ci[4, ] <- unlist(results_classic[[4]][1, 1:2])

pm10_sy <- matrix(data = NA, nrow = 4, ncol = 2)
pm10_sy[1, ] <- unlist(results_synth[[1]][1, 1:2])
pm10_sy[2, ] <- unlist(results_synth[[2]][1, 1:2])
pm10_sy[3, ] <- unlist(results_synth[[3]][1, 1:2])
pm10_sy[4, ] <- unlist(results_synth[[4]][1, 1:2])

library("gplots")
pdf(file = "chicago_synth.pdf", width = 6, height = 4)
par(mar = c(4,4,4,1))
plotCI(x = 1:4, y = pm10_ci[, 1], ui = pm10_ci[, 1] + 2 * pm10_ci[, 2], 
                                   li = pm10_ci[, 1] - 2 * pm10_ci[, 2],
       ylim = c(-1e-4, 1.5e-3), xlim = c(0.8, 4.5), 
       lwd = 2, xaxt = 'n',
       xlab = "Lag of pm10tmean", ylab = "Death ~ PM10",
       main = "Chicago, Ill. 1987-2000")
abline(h = 0)
plotCI(x = 1:4, y = pm10_sy[, 1], ui = pm10_sy[, 1] + 2 * pm10_sy[, 2],
                                   li = pm10_sy[, 1] - 2 * pm10_sy[, 2],
       ylim = c(-1e-4, 1.5e-3), add = TRUE, col = "blue", lwd = 2, xaxt = 'n')
segments(x0 = 1.15, x1 = 1.15, y0 = pm10_ci[1, 1], y1 = pm10_sy[1, 1], lty = 2)
segments(x0 = 2.15, x1 = 2.15, y0 = pm10_ci[2, 1], y1 = pm10_sy[2, 1], lty = 2)
segments(x0 = 3.15, x1 = 3.15, y0 = pm10_ci[3, 1], y1 = pm10_sy[3, 1], lty = 2)
segments(x0 = 4.15, x1 = 4.15, y0 = pm10_ci[4, 1], y1 = pm10_sy[4, 1], lty = 2)
text(x = 1.4, y = (pm10_sy[1, 1] + pm10_ci[1, 1]) / 2, label = paste0(round(10000 * (pm10_sy[1, 1] - pm10_ci[1, 1]), 1), "e-4"))
text(x = 2.4, y = (pm10_sy[2, 1] + pm10_ci[2, 1]) / 2, label = paste0(round(10000 * (pm10_sy[2, 1] - pm10_ci[2, 1]), 1), "e-4"))
text(x = 3.4, y = (pm10_sy[3, 1] + pm10_ci[3, 1]) / 2, label = paste0(round(10000 * (pm10_sy[3, 1] - pm10_ci[3, 1]), 1), "e-4"))
text(x = 4.4, y = (pm10_sy[4, 1] + pm10_ci[4, 1]) / 2, label = paste0(round(10000 * (pm10_sy[4, 1] - pm10_ci[4, 1]), 1), "e-4"))
axis(side = 1, line = 0, at = 1:4, labels = c("Lag-0", "Lag-1", "Lag-2", "Lag-3"))
dev.off()

#
#  Do it in ggplot2
#
library("ggplot2")
dat <- data.frame(Estimate = c(pm10_ci[, 1], pm10_sy[, 1]),
                  SE = c(pm10_ci[, 2], pm10_sy[, 2]), 
                  lagnum = 1:4,
                  ymin = c(pm10_ci[, 1] - 2 * pm10_ci[, 2],
                           pm10_sy[, 1] - 2 * pm10_sy[, 2]),
                  ymax = c(pm10_ci[, 1] + 2 * pm10_ci[, 2],
                           pm10_sy[, 1] + 2 * pm10_sy[, 2]), 
                  mod = c(rep("Sng", 4), rep("Syn", 4)))
lab_diff <- c(paste0(round(10000 * (dat$Estimate[5] - dat$Estimate[1]), 1), "e-4"),
              paste0(round(10000 * (dat$Estimate[6] - dat$Estimate[2]), 1), "e-4"),
              paste0(round(10000 * (dat$Estimate[7] - dat$Estimate[3]), 1), "e-4"),
              paste0(round(10000 * (dat$Estimate[8] - dat$Estimate[4]), 1), "e-4"))

dfs <- 8
epsilon <- 1e-7
bf.epsilon <- 1e-7
max_iter <- 100
bf.max_iter <- 100
resp <- "death"
airpoll <- "pm10tmean + l1pm10tmean + l2pm10tmean + l3pm10tmean"
tempVar <- "tmpd"
N <- dim(data.spec)[1]
timeBasis_ns <- ns(data.spec[, "Time"], df = dfs * length(unique(data.spec[, "Year"]))) 
col_vect <- which(apply(timeBasis_ns, MAR = 2, FUN = max) > 1e-3)
timeBasis_ns <- timeBasis_ns[, col_vect]

gam.stmt <- c("gam(as.formula(", paste0("), data = data.spec, family = quasi(log, mu), 
                      na.action = na.exclude, epsilon = ", epsilon, ", bf.epsilon = ", bf.epsilon, ", maxit = ", max_iter,
                      ", bf.maxit = ", bf.max_iter, ")"))
model.stmt1 <- paste0(resp, " ~ ", airpoll, " + as.factor(dow) + ",
                          tempVar, " + timeBasis_ns")
model.stmt2 <- paste0(resp, " ~ ", airpoll, " + as.factor(dow) + ",
                          tempVar, " + timeBasis_slp")

fit <- eval(parse(text = paste0(gam.stmt[1], "model.stmt1", gam.stmt[2])))
total_coef <- sum(summary.glm(fit)$coefficients[2:5, "Estimate"])
ci <- ggplot(dat, aes(x = lagnum, y = Estimate, colour = mod)) + 
        geom_segment(aes(x = 0.5, xend = 5, y = total_coef, yend = total_coef, 
                     colour = "black", alpha = 0.5), linetype = 2) + 
        xlab("Lag of PM10") + ylab("Association: Death ~ PM10") + 
        geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1) + 
        geom_point() + ylim(c(-1e-4, 1.5e-3)) + 
        xlim(c(0.5, 4.5)) +
        annotate("text", label = lab_diff, x = seq(1.4, 4.4, 1), 
                 y = (dat$Estimate[1:4] + dat$Estimate[5:8]) / 2) + 
        geom_segment(aes(x = 1.15, xend = 1.15, y = dat$Estimate[1], 
                     yend = dat$Estimate[5], colour = "black"), linetype = 1) + 
        geom_segment(aes(x = 2.15, xend = 2.15, y = dat$Estimate[2], 
                     yend = dat$Estimate[6], colour = "black"), linetype = 1) + 
        geom_segment(aes(x = 3.15, xend = 3.15, y = dat$Estimate[3], 
                     yend = dat$Estimate[7], colour = "black"), linetype = 1) + 
        geom_segment(aes(x = 4.15, xend = 4.15, y = dat$Estimate[4], 
                     yend = dat$Estimate[8], colour = "black"), linetype = 1) +
        theme(legend.position = "none")

pdf(file = "chicago_synth_compare.pdf", width = 9, height = 5)
par(mar = c(4,4,4,4))
plot(ci)
dev.off()

