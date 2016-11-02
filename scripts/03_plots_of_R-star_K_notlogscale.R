library(lubridate)

#### Step 1: Read in data ####
Rstar_Augexpt<-read.csv("data-processed/cellcount_Augexpt.csv") #16 and 25 degree experiments
Rstar_Septexpt<-read.csv("data-processed/cellcount_Septexpt.csv") #16 and 25 degree experiments


Rstar_int<-read.csv("data-processed/cell_count_init.csv") #dummy initial time data

head(Rstar_Augexpt)
str(Rstar_int)

#### Step 2: make date and time into readable time #### 
str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
Rstar_Septexpt$start_time<-ymd_hms(Rstar_Septexpt$start_time)
Rstar_int$start_time<-ymd_hm (Rstar_int$start_time)


#### Step 3: merge datasets and dummy initial data #### 
Rstar_expt<-rbind(Rstar_Augexpt, Rstar_Septexpt) #leave out initial as I want to plot biovolume
#Rstar_expt<-rbind(Rstar_Augexpt, Rstar_Septexpt, Rstar_int)

#### Step 4: calculate biovolume ####
Rstar_expt$biovolume<-Rstar_expt$cell_count*Rstar_expt$volume

Rstar_expt10<-subset(Rstar_expt, Rstar_expt$temperature==10)
Rstar_expt16<-subset(Rstar_expt, Rstar_expt$temperature==16)
Rstar_expt25<-subset(Rstar_expt, Rstar_expt$temperature==25)
Rstar_expt30<-subset(Rstar_expt, Rstar_expt$temperature==30)

Rstar_expt10$start_time
#### Step 3: plot growth and nitrate for both species by temperature####
pdf(file="./figures/Growth_nitrate_biovolume_across_time_and_temp_notlog.pdf", width = 8, height = 8)
par(mfrow=c(3,4))
par(mar=c(3,3,0.5, 0.5), oma=c(0,0,0,0))
with(Rstar_expt10, plot((cell_count)~start_time, type="n", main="10", ylab="", xlab="", ylim=c(0, 300000), xlim=c(min(start_time), max(start_time))))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points((cell_count)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points((cell_count)~start_time, col=2))
mtext(side = 2, line = 2, 'cells/ml')
legend("topleft", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot((cell_count)~start_time, ylab="", xlab="", type="n", main="16", ylim=c(0, 300000)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points((cell_count)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points((cell_count)~start_time, col=2))

with(Rstar_expt25, plot((cell_count)~start_time, ylab="", xlab="", type="n", main="25", ylim=c(0, 300000)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points((cell_count)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points((cell_count)~start_time, col=2))

with(Rstar_expt30, plot((cell_count)~start_time, ylab="", xlab="", type="n", main="30", ylim=c(0, 300000)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points((cell_count)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points((cell_count)~start_time, col=2))

(Rstar_expt25$start_time)
#



with(Rstar_expt10, plot((volume)~start_time, type="n", main="10", ylab="", xlab="", ylim=c(0, 1500)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt10, Rstar_expt10$species=="BB"), points((volume)~start_time, col=3))
mtext(side = 2, line = 2, 'cell size')
legend("topright", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot((volume)~start_time, ylab="", xlab="", type="n", main="16", ylim=c(0, 1500)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt16, Rstar_expt16$species=="BB"), points((volume)~start_time, col=3))


with(Rstar_expt25, plot((volume)~start_time, ylab="", xlab="", type="n", main="25", ylim=c(0, 1500)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt25, Rstar_expt25$species=="BB"), points((volume)~start_time, col=3))

with(Rstar_expt30, plot((volume)~start_time, ylab="", xlab="", type="n", main="30", ylim=c(0, 1500)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt30, Rstar_expt30$species=="BB"), points((volume)~start_time, col=3))


#
# biovolume data

with(Rstar_expt10, plot((biovolume)~start_time, type="n", ylab="", xlab="", ylim=c(0, 1.5E8)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points((biovolume)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points((biovolume)~start_time, col=2))
mtext(side = 2, line = 2, 'biovolume')
legend("topleft", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot((biovolume)~start_time, ylab="", xlab="", type="n", ylim=c(0, 1.5E8)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points((biovolume)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points((biovolume)~start_time, col=2))

with(Rstar_expt25, plot((biovolume)~start_time, ylab="", xlab="", type="n", ylim=c(0, 1.5E8)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points((biovolume)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points((biovolume)~start_time, col=2))

with(Rstar_expt30, plot((biovolume)~start_time, ylab="", xlab="", type="n", ylim=c(0, 1.5E8)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points((biovolume)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points((biovolume)~start_time, col=2))


sept13 - OCt28
aug12-Oct4
#add nitrate data
#
#par(mfrow=c(1,4))
#par(mar=c(3,3,3, 0.5), oma=c(0,0,0,0))

with(nitratedata10, plot(abs~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 0.2)))
with(subset(nitratedata10, nitratedata10$species=="CH"), points(abs~date_N_assay, col=1, pch=2))
with(subset(nitratedata10, nitratedata10$species=="TT"), points(abs~date_N_assay, col=2, pch=2))
mtext(side = 2, line = 2, 'Nitrate, uM')
mtext(side = 3, line = 1.5, '10°C')
legend("topright", pch=1, col=c(1,2), c("CH", "TT"))

with(nitratedata16, plot(abs~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 0.2)))
with(subset(nitratedata16, nitratedata16$species=="CH"), points(abs~date_N_assay, col=1, pch=2))
with(subset(nitratedata16, nitratedata16$species=="TT"), points(abs~date_N_assay, col=2, pch=2))
mtext(side = 3, line = 1.5, '16°C')

with(nitratedata25, plot(abs~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 0.2)))
with(subset(nitratedata25, nitratedata25$species=="CH"), points(abs~date_N_assay, col=1, pch=2))
with(subset(nitratedata25, nitratedata25$species=="TT"), points(abs~date_N_assay, col=2, pch=2))
mtext(side = 3, line = 1.5, '25°C')

with(nitratedata30, plot(abs~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 0.2)))
with(subset(nitratedata30, nitratedata30$species=="CH"), points(abs~date_N_assay, col=1, pch=2))
with(subset(nitratedata30, nitratedata30$species=="TT"), points(abs~date_N_assay, col=2, pch=2))
mtext(side = 3, line = 1.5, '30°C')


dev.off()

subset(Rstar_expt30, Rstar_expt30$volume==max(Rstar_expt30$volume))
#
# size data
#
par(mfrow=c(1,4))
par(mar=c(3,3,0.5, 0.5), oma=c(0,0,0,0))
with(Rstar_expt10, plot((volume)~start_time, type="n", main="10", ylab="", xlab="", ylim=c(0, 1500)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt10, Rstar_expt10$species=="BB"), points((volume)~start_time, col=3))
mtext(side = 2, line = 2, 'cell size')
legend("topright", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot((volume)~start_time, ylab="", xlab="", type="n", main="16", ylim=c(0, 1500)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt16, Rstar_expt16$species=="BB"), points((volume)~start_time, col=3))


with(Rstar_expt25, plot((volume)~start_time, ylab="", xlab="", type="n", main="25", ylim=c(0, 1500)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt25, Rstar_expt25$species=="BB"), points((volume)~start_time, col=3))

with(Rstar_expt30, plot((volume)~start_time, ylab="", xlab="", type="n", main="30", ylim=c(0, 1500)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points((volume)~start_time, col=2))
with(subset(Rstar_expt30, Rstar_expt30$species=="BB"), points((volume)~start_time, col=3))

#### Step 3: plot growth AND nitrate for both species by temperature####
par(mfrow=c(1,4))
with(Rstar_expt10, plot(log(cell_count)~start_time, type="n", main="10", ylab="", xlab="", ylim=c(0, 14)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points(log(cell_count)~start_time, col=2))
mtext(side = 2, line = 2, 'log cells/ml')
legend("bottomleft", pch=1, col=c(1,2), c("CH", "TT"))
mtext("Time", 1,3)

par(new = T)
with(nitratedata10, plot(abs~date_N_assay, type="n", axes=F, ylab="", xlab="", ylim=c(0, 0.2)))
with(subset(nitratedata10, nitratedata10$species=="CH"), points(abs~date_N_assay, col=1, pch=2))
with(subset(nitratedata10, nitratedata10$species=="TT"), points(abs~date_N_assay, col=2, pch=2))
axis(side = 4)
mtext(side = 4, line = 2, 'Nitrate')
mtext("Time", 1,3)

with(Rstar_expt16, plot(log(cell_count)~start_time, ylab="", xlab="", type="n", main="16", ylim=c(1, 14)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points(log(cell_count)~start_time, col=2))
mtext(side = 2, line = 2, 'log cells/ml')

par(new = T)
with(nitratedata16, plot(abs~date_N_assay, type="n", ylab="", xlab="", axes=F, ylim=c(0, 0.2)))
with(subset(nitratedata16, nitratedata16$species=="CH"), points(abs~date_N_assay, col=1, pch=2))
with(subset(nitratedata16, nitratedata16$species=="TT"), points(abs~date_N_assay, col=2, pch=2))
axis(side = 4)
mtext(side = 4, line = 3, 'Nitrate')
mtext("Time", 1,3)

with(Rstar_expt25, plot(log(cell_count)~start_time, ylab="", xlab="", type="n", main="25", ylim=c(1, 14)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points(log(cell_count)~start_time, col=2))
mtext(side = 2, line = 2, 'log cells/ml')

par(new = T)
with(nitratedata25, plot(abs~date_N_assay, type="n", ylab="", xlab="", axes=F, ylim=c(0, 0.2)))
with(subset(nitratedata25, nitratedata25$species=="CH"), points(abs~date_N_assay, col=1))
with(subset(nitratedata25, nitratedata25$species=="TT"), points(abs~date_N_assay, col=2))
axis(side = 4)
mtext(side = 4, line = 3, 'Nitrate')
mtext("Time", 1,3)

with(Rstar_expt30, plot(log(cell_count)~start_time, ylab="", xlab="", type="n", main="30", ylim=c(1, 14)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points(log(cell_count)~start_time, col=2))
mtext(side = 2, line = 2, 'log cells/ml')

par(new = T)
with(nitratedata30, plot(abs~date_N_assay, type="n", ylab="", xlab="", main="30", ylim=c(0, 0.2)))
with(subset(nitratedata30, nitratedata30$species=="CH"), points(abs~date_N_assay, col=1))
with(subset(nitratedata30, nitratedata30$species=="TT"), points(abs~date_N_assay, col=2))
axis(side = 4)
mtext(side = 4, line = 3, 'Nitrate')
mtext("Time", 1,3)



#R* for each species by temp
Rstar10<-subset(nitratedata10, nitratedata10$date_N_assay>"2016-09-30" & nitratedata10$species %in% c("TT", "CH"))
with(Rstar10, plot(nitrate~species))

Rstar16<-subset(nitratedata16, nitratedata16$date_N_assay>"2016-08-13" & nitratedata16$date_N_assay<"2016-09-22" & nitratedata16$species %in% c("TT", "CH"))
with(Rstar16, plot(nitrate~species))

Rstar25<-subset(nitratedata25, nitratedata25$date_N_assay>"2016-08-13" & nitratedata25$species %in% c("TT", "CH"))
with(Rstar25, plot(nitrate~species))

Rstar30<-subset(nitratedata30, nitratedata30$date_N_assay>"2016-08-15"& nitratedata30$species %in% c("TT", "CH"))
with(Rstar30, plot(nitrate~species))

par(mfrow=c(1,1))
par(mar=c(4,4,0.5, 0.5), oma=c(0,0,0,0))
Rstars<-rbind(Rstar10, Rstar16, Rstar25, Rstar30)
with(Rstars, plot(nitrate~temp, col=ifelse(species=="TT", 2, 1), xlab="", ylab=""))
legend("topright", pch=1, col=c(1,2), c("CH", "TT"))
mtext("R-star, uM Nitrate", 2, 3)
mtext("Temperature", 1, 3)

#learn from Joey how she plot means so nicely#
####

plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="CH")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))

plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="BB")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))
