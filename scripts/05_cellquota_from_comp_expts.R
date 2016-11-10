#cell quota is final N-init N / final cells - init cells
#do this for each temp, each species

library(lubridate)

#### read in data #
nitratedata <- read.csv("./data-processed/nitratesamples.csv") ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
Rstar_Augexpt<-read.csv("data-processed/cellcount_Augexpt.csv") #16 and 25 degree experiments
Rstar_Septexpt<-read.csv("data-processed/cellcount_Septexpt.csv") #16 and 25 degree experiments

#### Step 2: make date and time into readable time #### 
#str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
Rstar_Septexpt$start_time<-ymd_hms(Rstar_Septexpt$start_time)
nitratedata$date_N_assay<-ymd(nitratedata$date_N_assay)

#calculate medium N from a mean for N init
mediumnitrate<-subset(nitratedata, nitratedata$species=="Medium")
Ninit<-max(mediumnitrate$nitrate)

#for N final, use the FIRST N assay day, these are far enough away from 0 that they may be more accurate.
Nfinal10CH<-subset(nitratedata10, nitratedata10$date_N_assay=="2016-09-14"&species=="CH")
Nfinal10TT<-subset(nitratedata10, nitratedata10$date_N_assay=="2016-09-14"&species=="TT")
Nfinal30CH<-subset(nitratedata30, nitratedata30$date_N_assay=="2016-09-14"&species=="CH")
Nfinal30TT<-subset(nitratedata30, nitratedata30$date_N_assay=="2016-09-14"&species=="TT")
Nfinal16CH<-subset(nitratedata16, nitratedata16$date_N_assay=="2016-08-11"&species=="CH")
Nfinal16TT<-subset(nitratedata16, nitratedata16$date_N_assay=="2016-08-11"&species=="TT")
Nfinal25CH<-subset(nitratedata25, nitratedata25$date_N_assay=="2016-08-11"&species=="CH")
Nfinal25TT<-subset(nitratedata25, nitratedata25$date_N_assay=="2016-08-11"&species=="TT")

#init cells
cellsinit<-1000

#for cells final, use the final closest to the first N assay day, or interpolate between 2 dates, or fit a growth curve
Rstar_expt<-rbind(Rstar_Augexpt, Rstar_Septexpt) #leave out initial as I want to plot biovolume
Rstar_expt10<-subset(Rstar_expt, Rstar_expt$temperature==10)
Rstar_expt16<-subset(Rstar_expt, Rstar_expt$temperature==16)
Rstar_expt25<-subset(Rstar_expt, Rstar_expt$temperature==25)
Rstar_expt30<-subset(Rstar_expt, Rstar_expt$temperature==30)


Rstar_expt16<-Rstar_expt16[complete.cases(Rstar_expt16),]
cellsfinal10CH<-subset(Rstar_expt10, Rstar_expt10$start_time>"2016-09-14"& Rstar_expt10$start_time<"2016-09-17"&species=="CH")
cellsfinal30CH<-subset(Rstar_expt30, Rstar_expt30$start_time>"2016-09-14"& Rstar_expt30$start_time<"2016-09-17"&species=="CH")
cellsfinal16CH<-subset(Rstar_expt16, Rstar_expt16$start_time>"2016-08-11"& Rstar_expt16$start_time<"2016-08-13"&species=="CH")
cellsfinal25CH<-subset(Rstar_expt25, Rstar_expt25$start_time>"2016-08-11"& Rstar_expt25$start_time<"2016-08-13"&species=="CH")
cellsfinal10TT<-subset(Rstar_expt10, Rstar_expt10$start_time>"2016-09-14"& Rstar_expt10$start_time<"2016-09-17"&species=="TT")
cellsfinal30TT<-subset(Rstar_expt30, Rstar_expt30$start_time>"2016-09-14"& Rstar_expt30$start_time<"2016-09-17"&species=="TT")
cellsfinal16TT<-subset(Rstar_expt16, Rstar_expt16$start_time>"2016-08-11"& Rstar_expt16$start_time<"2016-08-13"&species=="TT")
cellsfinal25TT<-subset(Rstar_expt25, Rstar_expt25$start_time>"2016-08-11"& Rstar_expt25$start_time<"2016-08-13"&species=="TT")


cellsfinal16CH$replicate

subset(Nfinal16CH, Nfinal16CH$rep%in%cellsfinal16CH$replicate)$nitrate
cellsfinal10CH$quota<-(Ninit-subset(Nfinal10CH, Nfinal10CH$rep%in%cellsfinal10CH$replicate)$nitrate)/((cellsfinal10CH$cell_count-500)*1000) #umol/L*ml/cells*1L/1000ml
cellsfinal16CH$quota<-(Ninit-subset(Nfinal16CH, Nfinal16CH$rep%in%cellsfinal16CH$replicate)$nitrate)/((cellsfinal16CH$cell_count-500)*1000)
cellsfinal25CH$quota<-(Ninit-subset(Nfinal25CH, Nfinal25CH$rep%in%cellsfinal25CH$replicate)$nitrate)/((cellsfinal25CH$cell_count-500)*1000)
cellsfinal30CH$quota<-(Ninit-subset(Nfinal30CH, Nfinal30CH$rep%in%cellsfinal30CH$replicate)$nitrate)/((cellsfinal30CH$cell_count-500)*1000)
cellsfinal10TT$quota<-(Ninit-subset(Nfinal10TT, Nfinal10TT$rep%in%cellsfinal10TT$replicate)$nitrate)/((cellsfinal10TT$cell_count-500)*1000)
cellsfinal16TT$quota<-(Ninit-subset(Nfinal16TT, Nfinal16TT$rep%in%cellsfinal16TT$replicate)$nitrate)/((cellsfinal16TT$cell_count-500)*1000)
cellsfinal25TT$quota<-(Ninit-subset(Nfinal25TT, Nfinal25TT$rep%in%cellsfinal25TT$replicate)$nitrate)/((cellsfinal25TT$cell_count-500)*1000)
cellsfinal30TT$quota<-(Ninit-subset(Nfinal30TT, Nfinal30TT$rep%in%cellsfinal30TT$replicate)$nitrate)/((cellsfinal30TT$cell_count-500)*1000)

#plot quota as a function of temp
par(mfrow=c(1,1))
par(mar=c(3,3,0.5, 0.5), oma=c(0,0,0,0))

with(equilib10, plot(1:4~temp, type="n", ylab="", xlab="", ylim=c(0, 0.1)))
with(cellsfinal10CH, points(quota~temperature, col=1))
with(cellsfinal16CH, points(quota~temperature, col=1))
with(cellsfinal25CH, points(quota~temperature, col=1))
with(cellsfinal30CH, points(quota~temperature, col=1))
with(cellsfinal10TT, points(quota~temperature, col=3))
with(cellsfinal16TT, points(quota~temperature, col=3))
with(cellsfinal25TT, points(quota~temperature, col=3))
with(cellsfinal30TT, points(quota~temperature, col=3))

Ninit<-50
#use these quota estimates to calculate R-star
CH10Kreps$calc.rstar<-Ninit-(CH10Kreps$cell_count*cellsfinal10CH$quota*1000) #mmol/cell*cell/ml*1000ml/L
CH16Kreps$calc.rstar<-Ninit-(CH16Kreps$cell_count*cellsfinal16CH$quota*1000) #mmol/cell*cell/ml*1000ml/L
CH25Kreps$calc.rstar<-Ninit-(CH25Kreps$cell_count*cellsfinal25CH$quota*1000) #mmol/cell*cell/ml*1000ml/L
CH30Kreps$calc.rstar<-Ninit-(CH30Kreps$cell_count*cellsfinal30CH$quota*1000) #mmol/cell*cell/ml*1000ml/L
TT10Kreps$calc.rstar<-Ninit-(TT10Kreps$cell_count*cellsfinal10TT$quota*1000) #mmol/cell*cell/ml*1000ml/L
TT16Kreps$calc.rstar<-Ninit-(TT16Kreps$cell_count*cellsfinal16TT$quota*1000) #mmol/cell*cell/ml*1000ml/L
TT25Kreps$calc.rstar<-Ninit-(TT25Kreps$cell_count*cellsfinal25TT$quota*1000) #mmol/cell*cell/ml*1000ml/L
TT30Kreps$calc.rstar<-Ninit-(TT30Kreps$cell_count*cellsfinal30TT$quota*1000) #mmol/cell*cell/ml*1000ml/L

with(equilib10, plot(1:4~temp, type="n", ylab="", xlab="", ylim=c(-1000, 1000)))
with(CH10Kreps, points(calc.rstar~temperature, col=1))
with(CH16Kreps, points(calc.rstar~temperature, col=1))
with(CH25Kreps, points(calc.rstar~temperature, col=1))
with(CH30Kreps, points(calc.rstar~temperature, col=1))
with(TT10Kreps, points(calc.rstar~temperature, col=3))
with(TT16Kreps, points(calc.rstar~temperature, col=3))
with(TT25Kreps, points(calc.rstar~temperature, col=3))
with(TT30Kreps, points(calc.rstar~temperature, col=3))
