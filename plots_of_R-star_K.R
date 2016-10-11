library(lubridate)

#### Step 1: Read in data ####
Rstar_Augexpt<-read.csv("data-processed/cellcount_Augexpt.csv") #16 and 25 degree experiments
head(Rstar_Augexpt)
str(Rstar_Augexpt)

Rstar_int<-read.csv("data-processed/cell_count_init.csv") #dummy initial time data
head(Rstar_int)
str(Rstar_int)


min(Rstar_Augexpt$start_time)
#### Step 2: make date and time into readable time #### 
str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
Rstar_int$start_time<-ymd_hms(Rstar_int$start_time)

Rstar_Augexpt<-rbind(Rstar_Augexpt, Rstar_int)


#### Step 3: plot it ####
min(Rstar_Augexpt$start_time)
par(mfrow=c(1,3))
plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, type="n", main="TT", xlim=c())
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="TT"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="TT"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("16°C", "25°C"))


plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="CH")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))

plot(log(Rstar_Augexpt$cell_count)~Rstar_Augexpt$start_time, ty="n", main="BB")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==16), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==25), points(log(cell_count)~start_time, col=2))
