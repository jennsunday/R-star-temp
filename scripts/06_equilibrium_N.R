library(lubridate)

#### read in data #
nitratedata <- read.csv("data-processed/nitratesamples.csv") ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
Rstar_Augexpt<-read.csv("data-processed/cellcount_Augexpt.csv") #16 and 25 degree experiments
Rstar_Septexpt<-read.csv("data-processed/cellcount_Septexpt.csv") #16 and 25 degree experiments


#### Step 2: make date and time into readable time #### 
#str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
Rstar_Septexpt$start_time<-ymd_hms(Rstar_Septexpt$start_time)
nitratedata$date_N_assay<-ymd(nitratedata$date_N_assay)
names(nitratedata)

nitrate25<-subset(nitratedata, temp=="25")
#### subset to data after K is clearly reached ####
finalN10CH<-subset(nitratedata, nitratedata$date_N_assay>"2016-10-12" & temp=="10" & species=="CH")
finalN16CH<-subset(nitratedata, nitratedata$date_N_assay>"2016-09-21" & temp=="16" & species=="CH")
finalN25CH<-subset(nitratedata, nitratedata$date_N_assay>"2016-09-21" & temp=="25" & species=="CH")
finalN30CH<-subset(nitratedata, nitratedata$date_N_assay>"2016-10-12" & temp=="30" & species=="CH")
finalN10TT<-subset(nitratedata, nitratedata$date_N_assay>"2016-10-12" & temp=="10" & species=="TT")
finalN16TT<-subset(nitratedata, nitratedata$date_N_assay>"2016-09-21" & temp=="16" & species=="TT")
finalN25TT<-subset(nitratedata, nitratedata$date_N_assay>"2016-09-21" & temp=="25" & species=="TT")
finalN30TT<-subset(nitratedata, nitratedata$date_N_assay>"2016-10-12" & temp=="30" & species=="TT")

#get means
CHmeanRstar<-c(mean(finalN10CH$nitrate), mean(finalN16CH$nitrate), mean(finalN25CH$nitrate), mean(finalN30CH$nitrate))
TTmeanRstar<-c(mean(finalN10TT$nitrate), mean(finalN16TT$nitrate), mean(finalN25TT$nitrate), mean(finalN30TT$nitrate))
temp=c(10, 16, 25, 30)
CHmedianRstar<-c(median(finalN10CH$nitrate), median(finalN16CH$nitrate), median(finalN25CH$nitrate), median(finalN30CH$nitrate))
TTmedianRstar<-c(median(finalN10TT$nitrate), median(finalN16TT$nitrate), median(finalN25TT$nitrate), median(finalN30TT$nitrate))


par(mfrow=c(1,1))
plot(finalN10CH$nitrate~finalN10CH$temp, xlim=c(8,32), ylim=c(-5, 40), ylab="", las=1)
points(finalN16CH$nitrate~finalN16CH$temp)
points(finalN25CH$nitrate~finalN25CH$temp)
points(finalN30CH$nitrate~finalN30CH$temp)
points(finalN10TT$nitrate~finalN10TT$temp, col=3)
points(finalN16TT$nitrate~finalN16TT$temp, col=3)
points(finalN25TT$nitrate~finalN25TT$temp, col=3)
points(finalN30TT$nitrate~finalN30TT$temp, col=3)
mtext('R*, mM', 2, 3)

#plot means
points(CHmedianRstar~temp, cex=2)
lines(CHmedianRstar~temp)
points(TTmedianRstar~temp, cex=2, col=3)
lines(TTmedianRstar~temp, col=3)

