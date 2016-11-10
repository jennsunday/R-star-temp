library(lubridate)

#### read in data #
nitratedata <- read.csv("./data-processed/nitratesamples.csv") ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
Rstar_Augexpt<-read.csv("data-processed/cellcount_Augexpt.csv") #16 and 25 degree experiments
Rstar_Septexpt<-read.csv("data-processed/cellcount_Septexpt.csv") #16 and 25 degree experiments

# quota at 16 from July 31
CHQ16 <- 2.2992e-07
TTQ16 <- 3.381176e-07

#### Step 2: make date and time into readable time #### 
#str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
Rstar_Septexpt$start_time<-ymd_hms(Rstar_Septexpt$start_time)
nitratedata$date_N_assay<-ymd(nitratedata$date_N_assay)

#### subset to data after K is clearly reached ####
equilib10<-subset(Rstar_expt10, Rstar_expt10$start_time>"2016-10-04")
equilib16<-subset(Rstar_expt16, Rstar_expt16$start_time>"2016-08-16")
equilib25<-subset(Rstar_expt25, Rstar_expt25$start_time>"2016-08-16")
equilib30<-subset(Rstar_expt30, Rstar_expt30$start_time>"2016-09-16")


# use mixed effects model to calculate mean K #
library(nlme)
CH10<-subset(equilib10, equilib10$species=="CH")
mod<-lme(cell_count~1, random=~1|replicate, data=CH10)
summary(mod)
meanCH10<-79879.84
sdCH10<-3248.017

CH16<-subset(equilib16, equilib16$species=="CH")
mod<-lme(cell_count~1, random=~1|replicate, data=CH16)
summary(mod)
meanCH16<-131712.1
sdCH16<-5959.314

CH25<-subset(equilib25, equilib25$species=="CH")
mod<-lme(cell_count~1, random=~1|replicate, data=CH25)
summary(mod)
meanCH25<-142300.2
sdCH25<-5868.267

CH30<-subset(equilib30, equilib30$species=="CH")
mod<-lme(cell_count~1, random=~1|replicate, data=CH30)
summary(mod)
meanCH30<-58732.34
sdCH30<-2817.993

TT10<-subset(equilib10, equilib10$species=="TT")
mod<-lme(cell_count~1, random=~1|replicate, data=TT10)
summary(mod)
meanTT10<-29860.88
sdTT10<-2372.111

TT16<-subset(equilib16, equilib16$species=="TT")
mod<-lme(cell_count~1, random=~1|replicate, data=TT16)
summary(mod)
meanTT16<-63509.52
sdTT16<-4262.371

TT25<-subset(equilib25, equilib25$species=="TT")
mod<-lme(cell_count~1, random=~1|replicate, data=TT25)
summary(mod)
meanTT25<-91491.76
sdTT25<-3069.728

TT30<-subset(equilib30, equilib30$species=="TT")
mod<-lme(cell_count~1, random=~1|replicate, data=TT30)
summary(mod)
meanTT30<-33823.14
sdTT30<-3078.561

CHmeanKs<-c(meanCH10, meanCH16, meanCH25, meanCH30)
CHsdKs<-c(sdCH10, sdCH16, sdCH25, sdCH30)
TTmeanKs<-c(meanTT10, meanTT16, meanTT25, meanTT30)
TTsdKs<-c(sdTT10, sdTT16, sdTT25, sdTT30)
temp<-c(10, 16, 25, 30)

# aggregate data from each rep
CH10Kreps<-aggregate(CH10, list(CH10$replicate), mean)
CH16Kreps<-aggregate(CH16, list(CH16$replicate), mean)
CH25Kreps<-aggregate(CH25, list(CH25$replicate), mean)
CH30Kreps<-aggregate(CH30, list(CH30$replicate), mean)
TT10Kreps<-aggregate(TT10, list(TT10$replicate), mean)
TT16Kreps<-aggregate(TT16, list(TT16$replicate), mean)
TT25Kreps<-aggregate(TT25, list(TT25$replicate), mean)
TT30Kreps<-aggregate(TT30, list(TT30$replicate), mean)

# calc R-star using quota for each rep
CH10Kreps$calc.rstar<-28.7294-(CH10Kreps$cell_count*CHQ16*1000) #mmol/cell*cell/ml*1000ml/L
CH16Kreps$calc.rstar<-28.7294-(CH16Kreps$cell_count*CHQ16*1000)
CH25Kreps$calc.rstar<-28.7294-(CH25Kreps$cell_count*CHQ16*1000)
CH30Kreps$calc.rstar<-28.7294-(CH30Kreps$cell_count*CHQ16*1000)
TT10Kreps$calc.rstar<-28.7294-(TT10Kreps$cell_count*TTQ16*1000)
TT16Kreps$calc.rstar<-28.7294-(TT16Kreps$cell_count*TTQ16*1000)
TT25Kreps$calc.rstar<-28.7294-(TT25Kreps$cell_count*TTQ16*1000)
TT30Kreps$calc.rstar<-28.7294-(TT30Kreps$cell_count*TTQ16*1000)

CHmeanRs<-c(mean(CH10Kreps$calc.rstar), mean(CH16Kreps$calc.rstar), mean(CH25Kreps$calc.rstar), mean(CH30Kreps$calc.rstar))
TTmeanRs<-c(mean(TT10Kreps$calc.rstar), mean(TT16Kreps$calc.rstar), mean(TT25Kreps$calc.rstar), mean(TT30Kreps$calc.rstar))

#### plot K as a function of temp - all data ####
par(mfrow=c(1,1))
par(mar=c(3,3,0.5, 0.5), oma=c(0,0,0,0))

with(equilib10, plot(1:4~temp, type="n", ylab="", xlab="", ylim=c(0, 300000)))
with(subset(equilib10, equilib10$species=="CH"), points((cell_count)~temperature, col=1))
with(subset(equilib10, equilib10$species=="TT"), points((cell_count)~temperature, col=3))
with(subset(equilib16, equilib16$species=="CH"), points((cell_count)~temperature, col=1))
with(subset(equilib16, equilib16$species=="TT"), points((cell_count)~temperature, col=3))
with(subset(equilib25, equilib25$species=="CH"), points((cell_count)~temperature, col=1))
with(subset(equilib25, equilib25$species=="TT"), points((cell_count)~temperature, col=3))
with(subset(equilib30, equilib25$species=="CH"), points((cell_count)~temperature, col=1))
with(subset(equilib30, equilib25$species=="TT"), points((cell_count)~temperature, col=3))


#### plot K as a function of temp - means from reps ####
par(mfrow=c(1,1))
par(mar=c(3,4,0.5, 0.5), oma=c(0,0,0,0))
with(equilib10, plot(1:4~temp, type="n", ylab="", xlab="", ylim=c(0, 200000), las=1))
with(CH10Kreps, points((cell_count)~temperature, col=1))
with(CH16Kreps, points((cell_count)~temperature, col=1))
with(CH25Kreps, points((cell_count)~temperature, col=1))
with(CH30Kreps, points((cell_count)~temperature, col=1))
with(TT10Kreps, points((cell_count)~temperature, col=3))
with(TT16Kreps, points((cell_count)~temperature, col=3))
with(TT25Kreps, points((cell_count)~temperature, col=3))
with(TT30Kreps, points((cell_count)~temperature, col=3))

points(CHmeanKs~temp, pch=1, cex=2, col=1)
lines(CHmeanKs~temp)
points(TTmeanKs~temp, pch=1, cex=2, col=3)
lines(TTmeanKs~temp, col=3)

legend("topright", col=c(1,3), c("Chlamy", "T.tetra"), bty="n", pch=1)


#### plot K-biomass as a function of temp - means from reps ####
par(mfrow=c(1,1))
par(mar=c(3,4,0.5, 0.5), oma=c(0,0,0,0))
with(equilib10, plot(1:4~temp, type="n", ylab="", xlab="", ylim=c(0, 100000000), las=1))
with(CH10Kreps, points((biovolume)~temperature, col=1))
with(CH16Kreps, points((biovolume)~temperature, col=1))
with(CH25Kreps, points((biovolume)~temperature, col=1))
with(CH30Kreps, points((biovolume)~temperature, col=1))
with(TT10Kreps, points((biovolume)~temperature, col=3))
with(TT16Kreps, points((biovolume)~temperature, col=3))
with(TT25Kreps, points((biovolume)~temperature, col=3))
with(TT30Kreps, points((biovolume)~temperature, col=3))

points(CHmeanKs~temp, pch=1, cex=2, col=1)
lines(CHmeanKs~temp)
points(TTmeanKs~temp, pch=1, cex=2, col=3)
lines(TTmeanKs~temp, col=3)

legend("topright", col=c(1,3), c("T.tetra", "Chlamy"), bty="n", pch=1)


#### plot R-star calculated as a function of temp - means from reps ####
par(mfrow=c(1,1))
par(mar=c(3,4,0.5, 0.5), oma=c(0,0,0,0))
with(equilib10, plot(1:4~temp, type="n", ylab="", xlab="", ylim=c(-10, 30), las=1))
with(CH10Kreps, points((calc.rstar)~temperature, col=1))
with(CH16Kreps, points((calc.rstar)~temperature, col=1))
with(CH25Kreps, points((calc.rstar)~temperature, col=1))
with(CH30Kreps, points((calc.rstar)~temperature, col=1))
with(TT10Kreps, points((calc.rstar)~temperature, col=3))
with(TT16Kreps, points((calc.rstar)~temperature, col=3))
with(TT25Kreps, points((calc.rstar)~temperature, col=3))
with(TT30Kreps, points((calc.rstar)~temperature, col=3))


points(CHmeanRs~temp, pch=1, cex=2, col=1)
lines(CHmeanRs~temp)
points(TTmeanRs~temp, pch=1, cex=2, col=3)
lines(TTmeanRs~temp, col=3)

legend("topright", col=c(1,3), c("Chlamy", "T.tetra"), bty="n", pch=1)

