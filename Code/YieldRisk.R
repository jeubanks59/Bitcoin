# Script: YieldRisk.R
# Author: Joshua L. Eubanks
# Date 1 Aug 2018

# Installing packages
library('forecast')


#Bringing in the data

setwd("../Data")
DailySpot <- read.csv("Daily5yrSpot.csv", header = TRUE)
FutJun <- read.csv("DailyFutureJun18.csv", header = TRUE)
yield <- read.csv("T-Bill.csv", header = TRUE)

# Cleaning the data
DailySpot$logClose <- log(DailySpot$Close)
DailySpot$Date <- as.Date(DailySpot$Date, format="%m/%d/%Y")
DailySpot <- DailySpot[,-c(2:5)]

FutJun$logJun <- log(FutJun$Close)
FutJun$Date <- as.Date(FutJun$Date, format="%m/%d/%Y")
FutJun <- FutJun[,-c(2:7)]

yield$Date <- as.Date(yield$Date, format="%m/%d/%Y")

# Bringing csvs together
Stuff <- merge(DailySpot, FutJun, by=c("Date"))
Stuff <- merge(Stuff, yield, by=c("Date"))

# Finding T-t
Stuff$T <- as.numeric(as.Date("2018-06-29")-Stuff$Date, units="days")
cmeJun <- Stuff[Stuff$T > 0,]

# Finding delta
cmeJun$num <- cmeJun$logClose - cmeJun$logJun
cmeJun$cme <- cmeJun$num/cmeJun$T
cmeJun$convience <- cmeJun$cme + cmeJun$Ask.Price/100

setwd("../Figures")
pdf("SpotvFut.pdf")
plot(cmeJun$Date, exp(cmeJun$logClose), type = "l",xlab = "Date", ylab = "Price")
lines(cmeJun$Date, exp(cmeJun$logJun), col = "blue")

legend(as.Date("2018-03-01"), 16000, legend=c("Spot Price", "June Futures Contract"),
       col=c("black", "blue"), lty=1, cex=0.8)
dev.off()

# Ornstein Ulenbeck process
arone <- arima(x = diff(cmeJun$logClose), order=c(1,0,0))
summary(arone)

pdf("acfSpot.pdf")
Acf(arone$residuals, ylim = c(-0.5,0.5)
    , main = "")
dev.off()

pdf("pacfSpot.pdf")
Pacf(arone$residuals, ylim = c(-0.5,0.5), main = "")
dev.off()

# Regression of convenience yield
y <- diff(cmeJun$convience)
x <- cmeJun$convience[-1]
model <- lm(y~x)
summary(model)

pdf("acfCon.pdf")
Acf(model$residuals, ylim = c(-0.5,0.5)
    , main = "")
dev.off()

pdf("pacfCon.pdf")
Pacf(model$residuals, ylim = c(-0.5,0.5)
    , main = "")
dev.off()

# Simulating Lambda
nsims <- 2000
rows <- length(cmeJun$logClose)
cols <- nsims
times <- nsims* rows
nas <- rep(NA,times)
X <- matrix(nas,nrow=rows,ncol=cols)
rho <- cor(arone$residuals, model$residuals)
set.seed(123457)
for(i in 1:length(cmeJun$logClose)){
cmeJun$kalpha <- rnorm(length(cmeJun$logClose),0.0001090,0.0001688)
cmeJun$k <- rnorm(length(cmeJun$logClose),-0.0031982,0.0104181)
cmeJun$alpha <- cmeJun$kalpha/cmeJun$k
sigmad <-0.0002124
sigmadsq <- sigmad^2
sigmassq <- 0.00343
sigmas <- sqrt(sigmassq)
cmeJun$lambda <- (cmeJun$logClose - cmeJun$logJun - cmeJun$convience*(1-exp(-cmeJun$k*cmeJun$T))/cmeJun$k + (cmeJun$Ask.Price - cmeJun$alpha + (sigmadsq/(cmeJun$k)^2)-(sigmas*sigmad*rho)/cmeJun$k)*cmeJun$T + (sigmadsq/4)*(1-exp(-cmeJun$k*cmeJun$T*2))/(cmeJun$k)^(3) + (cmeJun$kalpha + sigmad*sigmas*rho - (sigmadsq/cmeJun$k))*(1-exp(-cmeJun$k*cmeJun$T))/(cmeJun$k)^(2))/(1-exp(-cmeJun$k*cmeJun$T)-(cmeJun$T/cmeJun$k))
for(j in 1:nsims){
X[i,j] <- cmeJun$lambda[i]
set.seed(j)
    }
}


lambda <- as.vector(X)

# t-Test of lambda
t.test(lambda)

# Plotting the ecdf of lambda
pdf("Lambda.pdf")
plot(ecdf(lambda),verticals = TRUE, do.points = FALSE, main = "")
dev.off()

# Histogram of lambda
pdf("hist.pdf")
hist(lambda, main = "", xlab = "")
dev.off()





