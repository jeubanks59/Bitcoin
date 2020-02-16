# Script: Volatility.R
# Author: Joshua L. Eubanks
# Date: 1 Aug 2018
setwd("../Data")
thirmin <- read.csv("30Min6MonthSpot.csv", header = TRUE)
thirmin$Date <- strptime(thirmin$Date, format="%m/%d/%Y  %H:%M")
thirmin$logClose <- log(thirmin$Close)

setwd("../Figures")
pdf("30MinTrend.pdf")
plot(thirmin$Date,thirmin$logClose, type = "l"
	, main = ""
	, xlab = "Date"
	, ylab = "Logarithm of Spot Price")
dev.off()

pdf("PercentChange.pdf")
plot(thirmin$Date[-length(thirmin$Date)],diff(thirmin$logClose), type = "l"
	, main = ""
	, xlab = "Date"
	, ylab = "% Change")
dev.off()
set.seed(123457)

# Weekly
t <- 0:336  
sig2 <- var(diff(thirmin$logClose))

x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
x <- c(0, cumsum(x))


nsim <- 10000
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 
    1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))

pdf("Brownian.pdf")
plot(t, X[1, ], ylim = c(-0.5,0.5), type = "l"
	, main = ""
	, xlab = "Time"
	, ylab = "% Change")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)
dev.off()


# 2 hours
t <- 0:4 
sig2 <- var(diff(thirmin$logClose))
 
x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
x <- c(0, cumsum(x))

X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 
    1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))

pdf("2hrBrownian.pdf")
plot(t, X[1, ], ylim = c(-0.05,0.05), type = "l"
	, main = ""
	, xlab = "Time"
	, ylab = "% Change")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)
dev.off()


# Daily
t <- 0:48  
sig2 <- var(diff(thirmin$logClose))
x <- rnorm(n = length(t) - 1, sd = sqrt(sig2))
x <- c(0, cumsum(x))

nsim <- 10000
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 
    1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))

pdf("DayBrownian.pdf")
plot(t, X[1, ], ylim = c(-0.15,0.15), type = "l"
	, main = ""
	, xlab = "Time"
	, ylab = "% Change")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)
dev.off()
