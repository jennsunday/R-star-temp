library(lubridate)

#### read in data #

nitratedata <- read.csv("./data-processed/nitratesamples.csv") ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

head(nitratedata)


#make date and time into readable time
str(nitratedata)
nitratedata$date_N_assay<-ymd(nitratedata$date_N_assay)
str(nitratedata)

with(nitratedata, plot(nitrate~date_N_assay, col=species))
head(nitratedata)

nitratedata10<-subset(nitratedata, nitratedata$temp==10)
nitratedata16<-subset(nitratedata, nitratedata$temp==16)
nitratedata25<-subset(nitratedata, nitratedata$temp==25)
nitratedata30<-subset(nitratedata, nitratedata$temp==30)

mediumnitrate<-subset(nitratedata, nitratedata$species=="Medium")
init.nitrate<-subset(nitratedata, nitratedata$species=="Medium" & nitratedata$date_N_assay=="2016-08-18")

mean(mediumnitrate$nitrate)
min(mediumnitrate$date_N_assay)-5
plot(mediumnitrate$abs~mediumnitrate$date_N_assay)
#### Step 3: plot nitrate for both species by temperature####
par(mfrow=c(1,4))
par(mar=c(3,3,1, 0.5), oma=c(0,1,0,0))
with(nitratedata16, plot(nitrate~date_N_assay, type="n", ylim=c(0, 40), las=1, xlim=c(min(nitratedata10$date_N_assay)-5, max(nitratedata10$date_N_assay)+5)))
with(subset(nitratedata10, nitratedata10$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata10, nitratedata10$species=="TT"), points(nitrate~date_N_assay, col=3))

with(nitratedata16, plot(nitrate~date_N_assay, type="n", ylim=c(0, 40), las=1, xlim=c(min(nitratedata16$date_N_assay)-5, max(nitratedata16$date_N_assay)+5)))
with(subset(nitratedata16, nitratedata16$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata16, nitratedata16$species=="TT"), points(nitrate~date_N_assay, col=3))

with(nitratedata16, plot(nitrate~date_N_assay, type="n", ylim=c(0, 40), las=1, xlim=c(min(nitratedata16$date_N_assay)-5, max(nitratedata16$date_N_assay)+5)))
with(subset(nitratedata25, nitratedata25$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata25, nitratedata25$species=="TT"), points(nitrate~date_N_assay, col=3))

with(nitratedata16, plot(nitrate~date_N_assay, type="n", ylim=c(0, 40), las=1, xlim=c(min(nitratedata10$date_N_assay)-5, max(nitratedata10$date_N_assay)+5)))
with(subset(nitratedata30, nitratedata30$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata30, nitratedata30$species=="TT"), points(nitrate~date_N_assay, col=3))

