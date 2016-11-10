library(lubridate)

#### Step 1: Read in data ####
Rstar_Augexpt<-read.csv("data-processed/cellcount_Augexpt.csv") #16 and 25 degree experiments
Rstar_Septexpt<-read.csv("data-processed/cellcount_Septexpt.csv") #16 and 25 degree experiments
#Rstar_int<-read.csv("data-processed/cell_count_init.csv") #dummy initial time data

head(Rstar_int)
str(Rstar_int)

#### Step 2: make date and time into readable time #### 
#str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
Rstar_Septexpt$start_time<-ymd_hms(Rstar_Septexpt$start_time)
#Rstar_int$start_time<-ymd_hm (Rstar_int$start_time)

head(Rstar_Septexpt)
#### Step 3: merge datasets and dummy initial data #### 
Rstar_expt<-rbind(Rstar_Augexpt, Rstar_Septexpt) #leave out initial as I want to plot biovolume
#Rstar_expt<-rbind(Rstar_Augexpt, Rstar_Septexpt, Rstar_int)

#### Step 4: calculate biovolume ####
Rstar_expt$biovolume<-Rstar_expt$cell_count*Rstar_expt$volume

Rstar_expt10<-subset(Rstar_expt, Rstar_expt$temperature==10)
Rstar_expt16<-subset(Rstar_expt, Rstar_expt$temperature==16)
Rstar_expt25<-subset(Rstar_expt, Rstar_expt$temperature==25)
Rstar_expt30<-subset(Rstar_expt, Rstar_expt$temperature==30)


#### Step 3: plot growth and nitrate for both species by temperature####
pdf(file="./figures/Growth_nitrate_biovolume_across_time_and_temp.pdf", width = 8, height = 8)
par(mfrow=c(1,4))
par(mar=c(3,3,1, 0.5), oma=c(0,1,0,0))
with(Rstar_expt10, plot(log(cell_count)~start_time, type="n", ylab="", xlab="", ylim=c(1, 14), las=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points(log(cell_count)~start_time, col=3))
mtext(side = 2, line = 2, 'log cells/ml')
mtext(side = 3, line = -2, '10°C')
#legend("bottomleft", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot(log(cell_count)~start_time, ylab="", xlab="", type="n", ylim=c(1, 14), las=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points(log(cell_count)~start_time, col=3))
mtext(side = 3, line = -2, '16°C')


with(Rstar_expt25, plot(log(cell_count)~start_time, ylab="", xlab="", type="n", ylim=c(1, 14), las=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points(log(cell_count)~start_time, col=3))
mtext(side = 3, line = -2, '25°C')

with(Rstar_expt30, plot(log(cell_count)~start_time, ylab="", xlab="", type="n", ylim=c(1, 14), las=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points(log(cell_count)~start_time, col=3))
mtext(side = 3, line = -2, '30°C')

#add nitrate data
#

with(nitratedata10, plot(nitrate~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 40)))
with(subset(nitratedata10, nitratedata10$species=="CH"), points(nitrate~date_N_assay, col=1, pch=2))
with(subset(nitratedata10, nitratedata10$species=="TT"), points(nitrate~date_N_assay, col=2, pch=2))
mtext(side = 2, line = 2, 'Nitrate')

with(nitratedata16, plot(nitrate~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 40)))
with(subset(nitratedata16, nitratedata16$species=="CH"), points(nitrate~date_N_assay, col=1, pch=2))
with(subset(nitratedata16, nitratedata16$species=="TT"), points(nitrate~date_N_assay, col=2, pch=2))

with(nitratedata25, plot(nitrate~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 40)))
with(subset(nitratedata25, nitratedata25$species=="CH"), points(nitrate~date_N_assay, col=1, pch=2))
with(subset(nitratedata25, nitratedata25$species=="TT"), points(nitrate~date_N_assay, col=2, pch=2))

with(nitratedata30, plot(nitrate~date_N_assay, type="n", ylab="", xlab="", ylim=c(0, 40)))
with(subset(nitratedata30, nitratedata30$species=="CH"), points(nitrate~date_N_assay, col=1, pch=2))
with(subset(nitratedata30, nitratedata30$species=="TT"), points(nitrate~date_N_assay, col=2, pch=2))

#
#add size data
with(Rstar_expt10, plot(log(biovolume)~start_time, type="n", ylab="", xlab="", ylim=c(8, 20)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points(log(biovolume)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points(log(biovolume)~start_time, col=2))
mtext(side = 2, line = 2, 'log biovolume')
legend("bottomleft", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot(log(biovolume)~start_time, ylab="", xlab="", type="n", ylim=c(8, 20)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points(log(biovolume)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points(log(biovolume)~start_time, col=2))

with(Rstar_expt25, plot(log(biovolume)~start_time, ylab="", xlab="", type="n", ylim=c(8, 20)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points(log(biovolume)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points(log(biovolume)~start_time, col=2))

with(Rstar_expt30, plot(log(biovolume)~start_time, ylab="", xlab="", type="n", ylim=c(8, 20)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points(log(biovolume)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points(log(biovolume)~start_time, col=2))

dev.off()


#
# size data
#
par(mfrow=c(1,4))
par(mar=c(3,3,0.5, 0.5), oma=c(0,0,0,0))
with(Rstar_expt10, plot((volume)~start_time, type="n", main="10", ylab="", xlab="", ylim=c(0, 1500)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points((volume)~start_time, col=2))
mtext(side = 2, line = 2, 'log cell size')
legend("topright", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot((volume)~start_time, ylab="", xlab="", type="n", main="16", ylim=c(0, 1500)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points((volume)~start_time, col=2))

with(Rstar_expt25, plot((volume)~start_time, ylab="", xlab="", type="n", main="25", ylim=c(0, 1500)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points((volume)~start_time, col=2))

with(Rstar_expt30, plot((volume)~start_time, ylab="", xlab="", type="n", main="30", ylim=c(0, 1500)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points((volume)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points((volume)~start_time, col=2))

#### Step 3: plot growth AND nitrate for both species by temperature####
par(mfrow=c(1,4))
with(Rstar_expt10, plot(log(cell_count)~start_time, type="n", main="10", ylab="", xlab="", ylim=c(0, 14)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points(log(cell_count)~start_time, col=2))
mtext(side = 2, line = 2, 'log cells/ml')
legend("bottomleft", pch=1, col=c(1,2), c("CH", "TT"))
mtext("Time", 1,3)

par(new = T)
with(nitratedata10, plot(nitrate~date_N_assay, type="n", axes=F, ylab="", xlab="", ylim=c(0, 40)))
with(subset(nitratedata10, nitratedata10$species=="CH"), points(nitrate~date_N_assay, col=1, pch=2))
with(subset(nitratedata10, nitratedata10$species=="TT"), points(nitrate~date_N_assay, col=3, pch=2))
axis(side = 4)
mtext(side = 4, line = 2, 'Nitrate')
mtext("Time", 1,3)

with(Rstar_expt16, plot(log(cell_count)~start_time, ylab="", xlab="", type="n", main="16", ylim=c(1, 14)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points(log(cell_count)~start_time, col=3))
mtext(side = 2, line = 2, 'log cells/ml')

par(new = T)
with(nitratedata16, plot(nitrate~date_N_assay, type="n", ylab="", xlab="", axes=F, ylim=c(0, 40)))
with(subset(nitratedata16, nitratedata16$species=="CH"), points(nitrate~date_N_assay, col=1, pch=2))
with(subset(nitratedata16, nitratedata16$species=="TT"), points(nitrate~date_N_assay, col=3, pch=2))
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



plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="CH")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))

plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="BB")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))
