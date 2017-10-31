#
#  Analyze the computed risks from 1-
#
library("ggplot2")
z <- qnorm(1-(1-0.95)/2)

load("res_NMMAPS_dlnm.RDa")
load("res_NMMAPS_synth.RDa")
res_NMMAPS_synth <- models
res_NMMAPS_synth <- res_NMMAPS_synth[, c("CD", "AP", "Estimate.quasi", "Std.Dev.quasi")]
res_NMMAPS_synth <- exp(res_NMMAPS_synth)

res_NMMAPS <- res_NMMAPS[row.names(res_NMMAPS) %in% row.names(res_NMMAPS_synth), ]
res_NMMAPS2 <- res_NMMAPS_synth

#
#  res_NMMAPS is the DLNM results (as listed)
#  res_NMMAPS2 is the synthetic lag results
#
names(res_NMMAPS) <- c("dlnm_low", "dlnm", "dlnm_high")
names(res_NMMAPS2) <- c("synth_low", "synth", "synth_high")

ggdat <- data.frame(x = row.names(res_NMMAPS),
                    low = c(res_NMMAPS$dlnm_low, res_NMMAPS2$synth_low),
                    high = c(res_NMMAPS$dlnm_high, res_NMMAPS2$synth_high),
                    y = c(res_NMMAPS$dlnm, res_NMMAPS2$synth),
                    Model = as.factor(c(rep("dlnm", 21), rep("Synth", 21))))

p <- ggplot(ggdat, aes(x = x, y = y, colour = Model)) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     geom_point(aes(x = x, y = y)) +
     geom_errorbar(aes(ymin = low, ymax = high), width = 0.2) +
     xlab("City, USA") +
     ylab("Relative Risk (Increase of 1 unit of Ozone)") +
     ylim(ymin = min(ggdat$low), max(ggdat$high)) 

png("NMMAPS_dlnm_vs_synth.png", width = 1000, height = 500)
par(mar = c(4, 4, 1, 1))
plot(p)
dev.off()


