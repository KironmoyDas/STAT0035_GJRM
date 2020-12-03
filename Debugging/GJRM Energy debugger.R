packages <- c("magic", "matrixcalc", "VGAM", "survey", "trust", "VineCopula",
              "graphics", "stats", "utils", "grDevices", "ggplot2",
              "matrixStats", "mnormt", "gamlss.dist", "Rmpfr", "scam",
              "survival", "psych", "copula", "distrEx", "numDeriv",
              "trustOptim", "evd", "ismev", "mgcv", "sp")
lapply(packages, require, character.only = TRUE)
library(GJRM)

data("energy", package = "MSwM")
energy$t <- 1:dim(energy)[1] # create time variable

eq.mu.1 <- Price ~ s(t, k = 10) + s(Oil, k = 10) + s(Coal, k = 10)
eq.mu.2 <- Demand ~ s(t, k = 10) + s(Oil, k = 10) + s(Gas, k = 10) + 
                    s(Coal, k = 10)
eq.sigma2.1 <- ~ s(t, k = 10)
eq.sigma2.2 <- ~ s(t, k = 10) + s(Oil, k = 10) + s(Gas, k = 10)
eq.theta <- ~ s(t, k = 10)
fl <- list(eq.mu.1, eq.mu.2, eq.sigma2.1, eq.sigma2.2, eq.theta)

eq.mu.1 <- Price ~ s(t, k = 60) + s(Oil) + s(Coal)
eq.mu.2 <- Demand ~ s(t, k = 60) + s(Oil) + s(Gas) + s(Coal)
eq.sigma2.1 <- ~ s(t, k = 60)
eq.sigma2.2 <- ~ s(t, k = 60) + s(Oil) + s(Gas)
eq.theta <- ~ s(t, k = 60)
fl <- list(eq.mu.1, eq.mu.2, eq.sigma2.1, eq.sigma2.2, eq.theta)
# Note: Gaussian copula model is the quickest (around 12 minutes)
# and also the one selected by AIC and BIC
outN <- gjrm(fl, margins = c("N", "GU"),
                  data = energy, gamlssfit = TRUE, Model="B")
# outF <- copulaReg(fl, margins = c("N", "GU"), BivD = "F",
#                   data = energy, gamlssfit = TRUE)
# outAMH <- copulaReg(fl, margins = c("N", "GU"), BivD = "AMH",
#                     data = energy, gamlssfit = TRUE)
# outFGM <- copulaReg(fl, margins = c("N", "GU"), BivD = "FGM",
#                     data = energy, gamlssfit = TRUE)
# outC0C90 <- copulaReg(fl, margins = c("N", "GU"), BivD = "C0C90",
#                       data = energy, gamlssfit = TRUE)
# outC0C270 <- copulaReg(fl, margins = c("N", "GU"), BivD = "C0C270",
#                        data = energy, gamlssfit = TRUE)
# outC180C90 <- copulaReg(fl, margins = c("N", "GU"), BivD = "C180C90",
#                         data = energy, gamlssfit = TRUE)
# outC180C270 <- copulaReg(fl, margins = c("N", "GU"), BivD = "C180C270",
#                          data = energy, gamlssfit = TRUE)
# the Joe and Gumbel families were also tried and not
# reported here to avoid repetitions
# fl <- list(eq.mu.1, eq.mu.2, eq.sigma2.1, eq.sigma2.2, ~ 1, eq.theta)
# outT <- copulaReg(fl, margins = c("N", "GU"), BivD = "T",
#                   data = energy, gamlssfit = TRUE)
conv.check(outN)
# conv.check(outF)
# conv.check(outAMH)
# conv.check(outFGM)
# conv.check(outT)
# conv.check(outC0C90)

# conv.check(outC0C270)
# conv.check(outC180C90)
# conv.check(outC180C270)
# AIC(outN, outF, outAMH, outFGM, outT, outC0C90, outC0C270,
#     outC180C90, outC180C270)
# BIC(outN, outF, outAMH, outFGM, outT, outC0C90, outC0C270,
#     outC180C90, outC180C270)
tpc <- "Histogram and Density of Residuals"
postR <- post.check(outN, main = tpc, main2 = tpc)
par(mfrow = c(2, 2), cex.main = 1.3)
acf(energy$Price, main = "Price")
acf(postR$qr1, main = "Quantile residuals of Price")
acf(energy$Demand, main = "Demand")
acf(postR$qr2, main = "Quantile residuals of Demand")
ss <- summary(outN, n.sim = 1000)
# n.sim is the no. of simulated coefficient vectors from the posterior
# distribution of the estimated model parameters, which are used to
# calculate intervals for non-linear functions of the model parameters
ss
# estimated tau and theta for each observation
outN$tau
outN$theta
# produce plot for theta
require(xts)
theta <- ts(outN$theta, frequency = 261, start = c(2002, 1))
CItheta1 <- ts(ss$CItheta[, 1], frequency = 261, start = c(2002, 1))
CItheta2 <- ts(ss$CItheta[, 2], frequency = 261, start = c(2002, 1))
par(cex.axis = 1.4, cex.lab = 1.6)
plot(theta, type = "l", ylim = c(-0.85, 0.97),
     ylab = expression(hat(theta)),
     xlab = "year", lwd = 2)
lines(CItheta1, lty = 2); lines(CItheta2, lty = 2)
# produce plots for first margin of copula model
split.screen(figs = c(2, 1))
st <- predict(outN, eq = 1, se.fit = TRUE, type = "iterms")
tfit <- ts(st$fit[, 1], frequency = 261, start = c(2002, 1))
sfit1 <- ts(st$fit[, 1] - 1.9665*st$se.fit[, 1],
            frequency = 261, start = c(2002, 1))
sfit2 <- ts(st$fit[, 1] + 1.9665*st$se.fit[, 1],
            frequency = 261, start = c(2002, 1))
plot(tfit, ylim = c(-3.9, 4.4), ylab = "s(time, 57.71)", xlab = "year")
lines(sfit1, lty = 2); lines(sfit2, lty = 2)
split.screen(figs = c(1, 2), screen = 2)
plot(outN, eq = 1, select = 2, scale = 0, seWithMean = TRUE)
screen(4)
plot(outN, eq = 1, select = 3, scale = 0, seWithMean = TRUE)
# produce plot for sigma21
sigma21 <- ts(outN$sigma21, frequency = 261, start = c(2002, 1))
CIsigma211 <- ts(ss$CIsig21[, 1], frequency = 261, start = c(2002, 1))
CIsigma212 <- ts(ss$CIsig21[, 2], frequency = 261, start = c(2002, 1))
par(cex.axis = 1.4, cex.lab = 1.5)
plot.ts(sigma21, type = "l", ylim = c(0.02, 50), log = "y",
        ylab = expression(hat(sigma)[1]^2),
        xlab = "year", lwd = 2)
lines(CIsigma211, lty = 2); lines(CIsigma212, lty = 2)          
