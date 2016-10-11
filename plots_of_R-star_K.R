library(lubridate)

#### Step 1: Read in data ####
Rstar_Augexpt<-read.csv("data-processed/cellcount_Augexpt.csv") #16 and 25 degree experiments
Rstar_Septexpt<-read.csv("data-processed/cellcount_Septexpt.csv") #16 and 25 degree experiments


Rstar_int<-read.csv("data-processed/cell_count_init.csv") #dummy initial time data

head(Rstar_int)
str(Rstar_int)

#### Step 2: make date and time into readable time #### 
str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
Rstar_Septexpt$start_time<-ymd_hms(Rstar_Septexpt$start_time)
Rstar_int$start_time<-ymd_hm (Rstar_int$start_time)

#### Step 3: merge datasets and dummy initial data #### 
Rstar_expt<-rbind(Rstar_Augexpt, Rstar_Septexpt, Rstar_int)

Rstar_expt10<-subset(Rstar_expt, Rstar_expt$temperature==10)
Rstar_expt16<-subset(Rstar_expt, Rstar_expt$temperature==16)
Rstar_expt25<-subset(Rstar_expt, Rstar_expt$temperature==25)
Rstar_expt30<-subset(Rstar_expt, Rstar_expt$temperature==30)

#### Step 3: plot growth for both species by temperature####
par(mfrow=c(1,4))
with(Rstar_expt10, plot(log(cell_count)~start_time, type="n", main="10", ylim=c(1, 14)))
with(subset(Rstar_expt10, Rstar_expt10$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt10, Rstar_expt10$species=="TT"), points(log(cell_count)~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt16, plot(log(cell_count)~start_time, type="n", main="16", ylim=c(1, 14)))
with(subset(Rstar_expt16, Rstar_expt16$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt16, Rstar_expt16$species=="TT"), points(log(cell_count)~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt25, plot(log(cell_count)~start_time, type="n", main="25", ylim=c(1, 14)))
with(subset(Rstar_expt25, Rstar_expt25$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt25, Rstar_expt25$species=="TT"), points(log(cell_count)~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("CH", "TT"))

with(Rstar_expt30, plot(log(cell_count)~start_time, type="n", main="30", ylim=c(1, 14)))
with(subset(Rstar_expt30, Rstar_expt30$species=="CH"), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_expt30, Rstar_expt30$species=="TT"), points(log(cell_count)~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("CH", "TT"))





plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="CH")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))

plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="BB")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))
