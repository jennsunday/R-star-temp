#goal: Plot nitrate concentration as a function of time
#goal: Identify equilibrial nitrate for each temperature and species
#goal: take into account temporal "drift" of nitrate absorbance method.

#read in packages
library(lubridate)
library(ggplot2)
library(grid)
library(gridExtra)
library(nlme)
library(broom)

#read in processed data
data<-read.csv("data-processed/nitratedatacombined_error_evap.csv")

#make time vector into lubridated time
data$date<-ymd(data$date.sampled)
data$date.measured<-ymd(data$date.measured)

#remove visual outliers (nitrate values abherrently high, contaminated samples, low volume N assays that were abherrent)
levels(data$error_code)
cdata<-subset(data, !data$error_code %in% c("contaminated", "vis_outlier", "remove_low_vol", "discard"))

#remove evaporation-relevant samples (all hot) after time period when dilutions were maintained
#cdata<-subset(cdata, cdata$temp<20|cdata$date<"2017-05-19")
#data<-subset(cdata, cdata$temp<33|cdata$date<"2017-05-15")
dim(cdata)

#take average of 2 absorbances
cdata$meanabs<-(cdata$absorbance.1+cdata$absorbance.2)/2
#remove NAs
cdata<-subset(cdata, complete.cases(cdata$meanabs))

#convert absorbance to nitrate
#for data before July 2nd, use standard curve from 280617; 
#for data from after July 2nd, use standard curve from 110717 
#(or, consider interpolating standard curve for every day between standard curve creation)

cdata$nitrate<-ifelse(cdata$date.measured<"2017-07-02", 
                     (cdata$meanabs-coef(st.2806)[1])/coef(st.2806)[2], 
                     ifelse(cdata$date.measured<"2017-07-13" & cdata$date.measured>"2017-07-02", 
                            (cdata$meanabs-coef(st.1207)[1])/coef(st.1207)[2], 
                            (cdata$meanabs-coef(st.1407)[1])/coef(st.1407)[2]))


#convert date into day of experiment
sp1<-subset(cdata, cdata$species==1)
sp2<-subset(cdata, cdata$species==2)
sp3<-subset(cdata, cdata$species==3)
sp4<-subset(cdata, cdata$species==4)

sp1$day<-sp1$date-min(sp1$date)
sp2$day<-sp2$date-min(sp2$date)
sp3$day<-sp3$date-min(sp3$date)
sp4$day<-sp4$date-min(sp4$date)

cdata<-rbind(sp1, sp2, sp3, sp4)
cdata$day<-as.numeric(cdata$day)
summary(cdata$day)

#plot data
par(mfrow=c(1,2))
with (cdata, plot(meanabs~date))
with (cdata, plot(nitrate~date))

#
#
#
#

#plot by temp
par(mfrow=c(2,3))
for (i in c(3, 10, 17, 24, 31, 38)){
  tempdata<-subset(cdata, cdata$temp==i)
  with(tempdata, plot(nitrate~date, ylim=c(0, 12), col=i))
}

#plot by temp and species along date
par(mfcol=c(6,4))
par(mar=c(3,3,0,0), oma=c(1,1,0.5, 0.5))
for(j in 1:4){
  for (i in c(3, 10, 17, 24, 31, 38)){
    tempdata<-subset(cdata, cdata$temp==i & cdata$species==j)
    with(tempdata, plot(nitrate~date, ylim=c(0, 12), col=i))
   
  }
}

#according to notes, remove highest temp of species 1 before May 15
sp1<-subset(cdata, cdata$species==1)
sp2<-subset(cdata, cdata$species==2)
sp3<-subset(cdata, cdata$species==3)
sp4<-subset(cdata, cdata$species==4)
sp1c<-subset(sp1, sp1$temp<35 | sp1$date<"2017-05-15")
cdata<-rbind(sp1c, sp2, sp3, sp4)

#make all the initial concentrations the same
day0<-subset(cdata, cdata$day=="0")
day0$n.fixed<-mean(day0$nitrate)
daynot0<-subset(cdata, cdata$day!="0")
daynot0$n.fixed<-daynot0$nitrate
cdata<-rbind(day0, daynot0)

write_csv(cdata, "data-processed/processed_nitrate_data.csv")
