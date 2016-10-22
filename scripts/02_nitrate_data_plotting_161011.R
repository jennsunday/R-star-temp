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

plot(mediumnitrate$abs~mediumnitrate$date_N_assay)
#### Step 3: plot nitrate for both species by temperature####
par(mfrow=c(1,4))
with(nitratedata10, plot(abs~date_N_assay, type="n", main="10", ylim=c(0, 40)))
with(subset(nitratedata10, nitratedata10$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata10, nitratedata10$species=="TT"), points(nitrate~date_N_assay, col=2))

legend("topright", pch=1, col=c(1,2), c("CH", "TT"))

with(nitratedata16, plot(abs~date_N_assay, type="n", main="16", ylim=c(0, 40)))
with(subset(nitratedata16, nitratedata16$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata16, nitratedata16$species=="TT"), points(nitrate~date_N_assay, col=2))

with(nitratedata25, plot(abs~date_N_assay, type="n", main="25", ylim=c(0, 40)))
with(subset(nitratedata25, nitratedata25$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata25, nitratedata25$species=="TT"), points(nitrate~date_N_assay, col=2))

with(nitratedata30, plot(abs~date_N_assay, type="n", main="30", ylim=c(0, 40)))
with(subset(nitratedata30, nitratedata30$species=="CH"), points(nitrate~date_N_assay, col=1))
with(subset(nitratedata30, nitratedata30$species=="TT"), points(nitrate~date_N_assay, col=2))
