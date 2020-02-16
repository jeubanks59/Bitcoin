# Script: Contango.R
# Author: Joshua L. Eubanks
# Date: 1 Aug 2018

# Removing scientific notation
options(scipen = 999)

install.packages("forecast",repos='http://cran.us.r-project.org'
)

setwd("../Data")
# Importing csvs
DailySpot <- read.csv("Daily5yrSpot.csv", header = TRUE)
FutJan <- read.csv("DailyFutureJan18.csv", header = TRUE)
FutFeb <- read.csv("DailyFutureFeb18.csv", header = TRUE)
FutMar <- read.csv("DailyFutureMar18.csv", header = TRUE)
FutApr <- read.csv("DailyFutureApr18.csv", header = TRUE)
FutMay <- read.csv("DailyFutureMay18.csv", header = TRUE)
FutJun <- read.csv("DailyFutureJun18.csv", header = TRUE)
FutJul <- read.csv("DailyFutureJul18.csv", header = TRUE)
FutAug <- read.csv("DailyFutureAug18.csv", header = TRUE)
FutSep <- read.csv("DailyFutureSep18.csv", header = TRUE)
JanCBOE <- read.csv("JanCBOE.csv", header = TRUE) 
FebCBOE <- read.csv("FebCBOE.csv", header = TRUE)
MarCBOE <- read.csv("MarCBOE.csv", header = TRUE)
AprCBOE <- read.csv("AprCBOE.csv", header = TRUE)
MayCBOE <- read.csv("MayCBOE.csv", header = TRUE)
JunCBOE <- read.csv("JunCBOE.csv", header = TRUE)
JulCBOE <- read.csv("JulCBOE.csv", header = TRUE)
AugCBOE <- read.csv("AugCBOE.csv", header = TRUE)
SepCBOE <- read.csv("SepCBOE.csv", header = TRUE)
yield <- read.csv("T-Bill.csv", header = TRUE)


setwd("../Figures")
#Cleaning csvs
DailySpot$logClose <- log(DailySpot$Close)
DailySpot$Date <- as.Date(DailySpot$Date, format="%m/%d/%Y")
DailySpot <- DailySpot[DailySpot$Date >= as.Date("2017-12-11"),]
DailySpot <- DailySpot[,-c(2:5)]

FutJan$logJan <- log(FutJan$Close)
FutJan$Date <- as.Date(FutJan$Date, format="%m/%d/%Y")
FutJan <- FutJan[,-c(2:7)]
FutFeb$logFeb <- log(FutFeb$Close)
FutFeb$Date <- as.Date(FutFeb$Date, format="%m/%d/%Y")
FutFeb <- FutFeb[,-c(2:7)]
FutMar$logMar <- log(FutMar$Close)
FutMar$Date <- as.Date(FutMar$Date, format="%m/%d/%Y")
FutMar <- FutMar[,-c(2:7)]
FutApr$logApr <- log(FutApr$Close)
FutApr$Date <- as.Date(FutApr$Date, format="%m/%d/%Y")
FutApr <- FutApr[,-c(2:7)]
FutMay$logMay <- log(FutMay$Close)
FutMay$Date <- as.Date(FutMay$Date, format="%m/%d/%Y")
FutMay <- FutMay[,-c(2:7)]
FutJun$logJun <- log(FutJun$Close)
FutJun$Date <- as.Date(FutJun$Date, format="%m/%d/%Y")
FutJun <- FutJun[,-c(2:7)]
FutJul$logJul <- log(FutJul$Close)
FutJul$Date <- as.Date(FutJul$Date, format="%m/%d/%Y")
FutJul <- FutJul[,-c(2:7)]
FutAug$logAug <- log(FutAug$Close)
FutAug$Date <- as.Date(FutAug$Date, format="%m/%d/%Y")
FutAug <- FutAug[,-c(2:7)]
FutSep$logSep <- log(FutSep$Close)
FutSep$Date <- as.Date(FutSep$Date, format="%m/%d/%Y")
FutSep <- FutSep[,-c(2:7)]

JanCBOE$Date <- as.Date(JanCBOE$Date, format="%m/%d/%Y")
JanCBOE$JanlogLast <- log(JanCBOE$Last.Price)
JanCBOE <- JanCBOE[,-c(2:4)]
FebCBOE$Date <- as.Date(FebCBOE$Date, format="%m/%d/%Y")
FebCBOE$FeblogLast <- log(FebCBOE$Last.Price)
FebCBOE <- FebCBOE[,-c(2:4)]
MarCBOE$Date <- as.Date(MarCBOE$Date, format="%m/%d/%Y")
MarCBOE$MarlogLast <- log(MarCBOE$Last.Price)
MarCBOE <- MarCBOE[,-c(2:4)]
AprCBOE$Date <- as.Date(AprCBOE$Date, format="%m/%d/%Y")
AprCBOE$AprlogLast <- log(AprCBOE$Last.Price)
AprCBOE <- AprCBOE[,-c(2:4)]
MayCBOE$Date <- as.Date(MayCBOE$Date, format="%m/%d/%Y")
MayCBOE$MaylogLast <- log(MayCBOE$Last.Price)
MayCBOE <- MayCBOE[,-c(2:4)]
JunCBOE$Date <- as.Date(JunCBOE$Date, format="%m/%d/%Y")
JunCBOE$JunlogLast <- log(JunCBOE$Last.Price)
JunCBOE <- JunCBOE[,-c(2:4)]
JulCBOE$Date <- as.Date(JulCBOE$Date, format="%m/%d/%Y")
JulCBOE$JullogLast <- log(JulCBOE$Last.Price)
JulCBOE <- JulCBOE[,-c(2:4)]
AugCBOE$Date <- as.Date(AugCBOE$Date, format="%m/%d/%Y")
AugCBOE$AuglogLast <- log(AugCBOE$Last.Price)
AugCBOE <- AugCBOE[,-c(2:4)]
SepCBOE$Date <- as.Date(SepCBOE$Date, format="%m/%d/%Y")
SepCBOE$SeplogLast <- log(SepCBOE$Last.Price)
SepCBOE <- SepCBOE[,-c(2:4)]

yield$Date <- as.Date(yield$Date, format="%m/%d/%Y")
yield <- yield[yield$Date >= as.Date("2017-12-11"),]


# Creating a big table
Stuff <- merge(DailySpot, FutJan, by=c("Date"), all.x=TRUE)
Stuff <- merge(Stuff, FutFeb, by=c("Date"), all.x=TRUE)
Stuff <- merge(Stuff, FutMar, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, FutApr, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, FutMay, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, FutJun, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, FutJul, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, FutAug, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, FutSep, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, JanCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, FebCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, MarCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, AprCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, MayCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, JunCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, JulCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, AugCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, SepCBOE, all.x=T, by=c("Date"))
Stuff <- merge(Stuff, yield, all.x=T, by=c("Date"))


# Making equation 1 for all months and exchanges

# Jan CME
Stuff$T <- as.numeric(as.Date("2018-01-26")-Stuff$Date, units="days")
cmeJan <- Stuff[Stuff$T > 0,]
cmeJan$num <- cmeJan$logClose - cmeJan$logJan
cmeJan$cme <- cmeJan$num/cmeJan$T
cmeJan <- na.omit(cmeJan$cme)

# Feb CME
Stuff$T <- as.numeric(as.Date("2018-02-23")-Stuff$Date, units="days")
cmeFeb <- Stuff[Stuff$T > 0,]
cmeFeb$num <- cmeFeb$logClose - cmeFeb$logFeb
cmeFeb$cme <- cmeFeb$num/cmeFeb$T
cmeFeb <- na.omit(cmeFeb$cme)

# Mar CME
Stuff$T <- as.numeric(as.Date("2018-03-29")-Stuff$Date, units="days")
cmeMar <- Stuff[Stuff$T > 0,]
cmeMar$num <- cmeMar$logClose - cmeMar$logMar
cmeMar$cme <- cmeMar$num/cmeMar$T
cmeMar <- na.omit(cmeMar$cme)

# Apr CME
Stuff$T <- as.numeric(as.Date("2018-04-27")-Stuff$Date, units="days")
cmeApr <- Stuff[Stuff$T > 0,]
cmeApr$num <- cmeApr$logClose - cmeApr$logApr
cmeApr$cme <- cmeApr$num/cmeApr$T
cmeApr <- na.omit(cmeApr$cme)

# May CME
Stuff$T <- as.numeric(as.Date("2018-05-25")-Stuff$Date, units="days")
cmeMay <- Stuff[Stuff$T > 0,]
cmeMay$num <- cmeMay$logClose - cmeMay$logMay
cmeMay$cme <- cmeMay$num/cmeMay$T
cmeMay <- na.omit(cmeMay$cme)

# Jun CME
Stuff$T <- as.numeric(as.Date("2018-06-29")-Stuff$Date, units="days")
cmeJun <- Stuff[Stuff$T > 0,]
cmeJun$num <- cmeJun$logClose - cmeJun$logJun
cmeJun$cme <- cmeJun$num/cmeJun$T
cmeJun <- na.omit(cmeJun$cme)

# Jul CME
Stuff$T <- as.numeric(as.Date("2018-07-27")-Stuff$Date, units="days")
cmeJul <- Stuff[Stuff$T > 0,]
cmeJul$num <- cmeJul$logClose - cmeJul$logJul
cmeJul$cme <- cmeJul$num/cmeJul$T
cmeJul <- na.omit(cmeJul$cme)

# Aug CME
Stuff$T <- as.numeric(as.Date("2018-08-31")-Stuff$Date, units="days")
cmeAug <- Stuff[Stuff$T > 0,]
cmeAug$num <- cmeAug$logClose - cmeAug$logAug
cmeAug$cme <- cmeAug$num/cmeAug$T
cmeAug <- na.omit(cmeAug$cme)

# Sep CME
Stuff$T <- as.numeric(as.Date("2018-08-28")-Stuff$Date, units="days")
cmeSep <- Stuff[Stuff$T > 0,]
cmeSep$num <- cmeSep$logClose - cmeSep$logSep
cmeSep$cme <- cmeSep$num/cmeSep$T
cmeSep <- na.omit(cmeSep$cme)


# Jan CBOE
Stuff$T <- as.numeric(as.Date("2018-01-17")-Stuff$Date, units="days")
cboeJan <- Stuff[Stuff$T > 0,]
cboeJan$num <- cboeJan$logClose - cboeJan$JanlogLast
cboeJan$cboe <- cboeJan$num/cboeJan$T
cboeJan <- na.omit(cboeJan$cboe)

# Feb cboe
Stuff$T <- as.numeric(as.Date("2018-02-14")-Stuff$Date, units="days")
cboeFeb <- Stuff[Stuff$T > 0,]
cboeFeb$num <- cboeFeb$logClose - cboeFeb$FeblogLast
cboeFeb$cboe <- cboeFeb$num/cboeFeb$T
cboeFeb <- na.omit(cboeFeb$cboe)

# Mar cboe
Stuff$T <- as.numeric(as.Date("2018-03-14")-Stuff$Date, units="days")
cboeMar <- Stuff[Stuff$T > 0,]
cboeMar$num <- cboeMar$logClose - cboeMar$MarlogLast
cboeMar$cboe <- cboeMar$num/cboeMar$T
cboeMar <- na.omit(cboeMar$cboe)

# Apr cboe
Stuff$T <- as.numeric(as.Date("2018-04-18")-Stuff$Date, units="days")
cboeApr <- Stuff[Stuff$T > 0,]
cboeApr$num <- cboeApr$logClose - cboeApr$AprlogLast
cboeApr$cboe <- cboeApr$num/cboeApr$T
cboeApr <- na.omit(cboeApr$cboe)

# May cboe
Stuff$T <- as.numeric(as.Date("2018-05-16")-Stuff$Date, units="days")
cboeMay <- Stuff[Stuff$T > 0,]
cboeMay$num <- cboeMay$logClose - cboeMay$MaylogLast
cboeMay$cboe <- cboeMay$num/cboeMay$T
cboeMay <- na.omit(cboeMay$cboe)

# Jun cboe
Stuff$T <- as.numeric(as.Date("2018-06-13")-Stuff$Date, units="days")
cboeJun <- Stuff[Stuff$T > 0,]
cboeJun$num <- cboeJun$logClose - cboeJun$JunlogLast
cboeJun$cboe <- cboeJun$num/cboeJun$T
cboeJun <- na.omit(cboeJun$cboe)

# Jul cboe
Stuff$T <- as.numeric(as.Date("2018-07-18")-Stuff$Date, units="days")
cboeJul <- Stuff[Stuff$T > 0,]
cboeJul$num <- cboeJul$logClose - cboeJul$JullogLast
cboeJul$cboe <- cboeJul$num/cboeJul$T
cboeJul <- na.omit(cboeJul$cboe)

# Aug cboe
Stuff$T <- as.numeric(as.Date("2018-08-15")-Stuff$Date, units="days")
cboeAug <- Stuff[Stuff$T > 0,]
cboeAug$num <- cboeAug$logClose - cboeAug$AuglogLast
cboeAug$cboe <- cboeAug$num/cboeAug$T
cboeAug <- na.omit(cboeAug$cboe)

# Sep cboe
Stuff$T <- as.numeric(as.Date("2018-09-19")-Stuff$Date, units="days")
cboeSep <- Stuff[Stuff$T > 0,]
cboeSep$num <- cboeSep$logClose - cboeSep$SeplogLast
cboeSep$cboe <- cboeSep$num/cboeSep$T
cboeSep <- na.omit(cboeSep$cboe)


# Making edf function
EDF <- function(y,vector){
indicator <- 0
    for( i in 1:length(vector)){
        if( y>= vector[i]){
            indicator <- indicator + 1
        }
        prob <- indicator/length(vector)
    }
    return(prob)
}


# Jan
inc <- 0.0000001
increments <- seq(-0.032
                , .006
                , inc)

CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeJan)
    CBOE[i] <- EDF(increments[i],cboeJan)
}



pdf("CMEvsCBOEperMonth.pdf")
par(mfrow=c(3,3))
plot(increments,CME, type = "l"
    , main = "January"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

# Feb
increments <- seq(-0.003
                , .0032
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeFeb)
    CBOE[i] <- EDF(increments[i],cboeFeb)
}
plot(increments,CME, type = "l"
    , main = "February"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

#Mar
increments <- seq(-0.003
                , .0042
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeMar)
    CBOE[i] <- EDF(increments[i],cboeMar)
}
plot(increments,CME, type = "l"
    , main = "March"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

#Apr
increments <- seq(-0.005
                , .03
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeApr)
    CBOE[i] <- EDF(increments[i],cboeApr)
}
plot(increments,CME, type = "l"
    , main = "April"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

#May
increments <- seq(-0.005
                , .004
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeMay)
    CBOE[i] <- EDF(increments[i],cboeMay)
}
plot(increments,CME, type = "l"
    , main = "May"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

#Jun
increments <- seq(-0.00125
                , .00125
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeJun)
    CBOE[i] <- EDF(increments[i],cboeJun)
}
plot(increments,CME, type = "l"
    , main = "June"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

#Jul
increments <- seq(-0.0038
                , .001
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeJul)
    CBOE[i] <- EDF(increments[i],cboeJul)
}
plot(increments,CME, type = "l"
    , main = "July"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

#Aug
increments <- seq(-0.006
                , .001
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeAug)
    CBOE[i] <- EDF(increments[i],cboeAug)
}
plot(increments,CME, type = "l"
    , main = "August"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")

#Sep
increments <- seq(-0.009
                , .001
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cmeSep)
    CBOE[i] <- EDF(increments[i],cboeSep)
}
plot(increments,CME, type = "l"
    , main = "September"
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")
dev.off()


#Overall

cme <- c(cmeJan,cmeFeb,cmeMar,cmeApr,cmeMay,cmeJul,cmeAug,cmeSep
,cmeJun)
cboe <- c(cboeJan,cboeFeb,cboeMar,cboeApr,cboeMay,cboeJun,cboeJul,cboeAug,cboeSep)

increments <- seq(-0.032
                , .01
                , inc)
CME  <- numeric(length(increments))
CBOE <- numeric(length(increments))
LowerCME <- numeric(length(increments))
UpperCME <- numeric(length(increments))
LowerCBOE <- numeric(length(increments))
UpperCBOE <- numeric(length(increments))
alpha <- 0.05

for(i in 1:length(increments)){
    CME[i] <- EDF(increments[i],cme)
    CBOE[i] <- EDF(increments[i],cboe)
    LowerCME[i] <- CME[i] - sqrt((log(2/alpha)/(2*length(CME))))
    if(LowerCME[i]<0){
        LowerCME[i] <- 0  
    }
    UpperCME[i] <- CME[i] + sqrt((log(2/alpha)/(2*length(CME))))
    if(UpperCME[i]>1){
        UpperCME[i] <- 1  
    }
    LowerCBOE[i] <- CBOE[i] - sqrt((log(2/alpha)/(2*length(CME))))
    if(LowerCBOE[i]<0){
        LowerCBOE[i] <- 0  
    }
    UpperCBOE[i] <-CBOE[i] + sqrt((log(2/alpha)/(2*length(CME))))
    if(UpperCBOE[i]>1){
        UpperCBOE[i] <- 1  
    }

}
pdf("Overall.pdf")
par(mfrow=c(1,1))
plot(increments,CME, type = "l"
    , main = ""
    , xlab = ""
    , ylab = "")
lines(increments,CBOE, col = "blue")
legend(-0.03, 0.8, legend=c("CME", "CBOE"),
       col=c("black", "blue"), lty=1, cex=0.8)
dev.off()
pdf("CME.pdf")
plot(increments,CME, type = "l"
    , main = ""
    , xlab = ""
    , ylab = "")
lines(increments,LowerCME, lty = "dashed", col = "blue")
lines(increments,UpperCME, lty = "dashed", col = "blue")
dev.off()
pdf("CBOE.pdf")
plot(increments,CBOE, type = "l"
    , main = ""
    , xlab = ""
    , ylab = "")
lines(increments,LowerCBOE, col = "blue", lty = "dashed")
lines(increments,UpperCBOE, col = "blue", lty = "dashed")
dev.off()

D <- max(abs(CME-CBOE))
n <- length(cme)
m <- length(cboe)
calc <- 1.36*sqrt((n+m)/(n*m))

ks.test(cme, cboe)

calc

t.test(cme, mu = 0)
t.test(cboe, mu = 0)