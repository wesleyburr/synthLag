#
#  Simulation 1
#  * independent WN(0, \sigma^2) base
#  * one periodic component
#  * equal \omega terms
#  * allow \phi^{(1)} to vary smoothly
#
set.seed(67)
N <- M <- 1000
tarr <- 1:1000
phi_var <- seq(-pi, pi, pi/100)
cors <- vector(length = N)
lmod <- vector(length = N)
res1 <- res2 <- vector(length = length(phi_var))

# Allow i to vary smoothly from -pi to pi in steps of pi/100
for(i in 1:length(phi_var)) {
    # Do 1000 runs for each offset
    for(j in 1:M) {
      wn1 <- rnorm(N, 0, 1)
      wn2 <- rnorm(N, 0, 1)
      x1 <- wn1 + 1.0 * cos(2 * pi * 1/4 * tarr + phi_var[i])
      x2 <- wn2 + 1.0 * cos(2 * pi * 1/4 * tarr + 0)
      lmod[j] <- lm(x2 ~ x1)$coefficients[2]
    }
    res1[i] <- mean(lmod)
    res2[i] <- sd(lmod)
}

library("ggplot2")
results <- data.frame(phi = phi_var, coeff = res1, ymin = res1 + qnorm(0.025) * res2,
                      ymax = res1 + qnorm(0.975) * res2)
p <- ggplot(results, aes(phi)) + 
       geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70") +
       geom_line(aes(y = coeff)) + 
       geom_hline(yintercept = c(-1/3, 0, 1/3), lty = c(2, 1, 2),
                  col = c("blue", "black", "blue")) + 
       ylab("Coefficient of Linear Model") +
       xlab(expression(phi))
plot(p)
pdf(file = "sim1.pdf", width = 8, height = 6)
par(mar = c(4, 4, 0.5, 0.5))
plot(p)
dev.off()

#
#  Simulation 2
#  * independent WN(0, \sigma^2) base
#  * two periodic components each
#  * equal \omega terms for each pairing
#  * allow \phi_1^{(1)} and \phi_2^{(1)} to vary smoothly
#
N <- 1000
M <- 1000
phi_var1 <- seq(-pi, pi, pi/20)
phi_var2 <- seq(-pi, pi, pi/20)
cors <- vector(length = N)
lmod <- vector(length = N)

comp <- function(X, Y) {
  for(j in 1:M) {
    wn1 <- rnorm(N, 0, 1)
    wn2 <- rnorm(N, 0, 1)
    x1 <- wn1 + cos(2 * pi * 1/4 * tarr + X) + cos(2 * pi * 1/6 * tarr + Y)
    x2 <- wn2 + cos(2 * pi * 1/4 * tarr + 0.0) + cos(2 * pi * 1/6 * tarr + 0.0)
    lmod[j] <- lm(x2 ~ x1)$coefficients[2]
  }
  c(mean(lmod), sd(lmod))
}
res1 <- res2 <- matrix(data = NA, nrow = length(phi_var1), ncol = length(phi_var2))

for(j in 1:length(phi_var1)) {
  cat(".")
  for(k in 1:length(phi_var2)) {
    vals <- comp(phi_var1[j], phi_var2[k])
    res1[j, k] <- vals[1]
    res2[j, k] <- vals[2]
  }
}
cat("\n")

library(reshape2)
colnames(res1) <- phi_var1
rownames(res1) <- phi_var2
melted_lmod <- melt(res1)

p <- ggplot(data = melted_lmod, aes(x=Var1, y=Var2, fill=value)) + 
     geom_tile(aes(fill = value), colour = "white") +
     scale_fill_gradient(low = "white", high = "steelblue") +
     ylab(expression(phi[2])) +
     xlab(expression(phi[1]))

pdf(file = "sim2.pdf", width = 8, height = 6)
par(mar = c(4, 4, 0.5, 0.5))
plot(p)
dev.off()


