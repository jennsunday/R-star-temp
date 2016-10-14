### using N standards from each date to determine N concentration #

library(lubridate)
library(dplyr)
nitratedata <- read.csv("./data/nitrate_data/r-star_Nitrate_data_161011.csv") ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

head(nitratedata)

standards<-subset(nitratedata, nitratedata$species=="Standard")[,c(1,2,5,8,9)]

standards

#make a standard curve for each data
standards$date_N_assay<-as.factor(standards$date_N_assay)
par(mfrow=c(3,4))
for (i in 1:length(levels(standards$date_N_assay))){
  stan<-subset(standards, standards$date_N_assay==levels(standards$date_N_assay)[i])
  with(stan, plot(abs~nitrate, main=date_N_assay[1], col=i, ylim=c(0, 0.15), xlim=c(0, 30)))
  mod<-lm(abs~nitrate, data=stan)
  abline(mod, col=i)
}

#See what occurs when plots overlap
par(mfrow=c(1,1))
standards$date_N_assay<-as.factor(standards$date_N_assay)
slopes<-rep("NA", length=length(levels(standards$date_N_assay)))
int<-rep("NA", length=length(levels(standards$date_N_assay)))
with(standards, plot(abs~nitrate, type="n"))
for (i in 1:length(levels(standards$date_N_assay))){
  stan<-subset(standards, standards$date_N_assay==levels(standards$date_N_assay)[i])
  with(stan, points(abs~nitrate, main=date_N_assay[1], col=i))
  mod<-lm(abs~nitrate, data=stan)
  #clip(min(stan$nitrate)*0.8, max(stan$nitrate)*1.1, min(stan$abs)*0.8, max(stan$abs)*1.1)
  abline(mod, col=i)
  slopes[i]<-coef(mod)[2]
  int[i]<-coef(mod)[1]
}


#plot slopes against date to see pattern in drift
#make date and time into readable time
str(standards)
standards$date_N_assay<-dmy(standards$date_N_assay)
str(standards)

plot(slopes~unique(standards$date_N_assay))
plot(int~unique(standards$date_N_assay))

#don't see a pattern, but I think the slope also changes with the range of Ns assayed.

#apply standard curve to absorbance data from each time point
standards$date_N_assay<-as.factor(standards$date_N_assay)#make date a factor again
nitratesamples<-subset(nitratedata, nitratedata$species!="Standard")
nitratesamples$date_N_assay<-dmy(nitratesamples$date_N_assay)
nitratesamples$date_N_assay<-as.factor(nitratesamples$date_N_assay)

for(i in 1:length(levels(standards$date_N_assay))){
  Nset<-subset(nitratesamples, nitratesamples$date_N_assay==levels(standards$date_N_assay)[i])
  stan<-subset(standards, standards$date_N_assay==levels(standards$date_N_assay)[i])
  mod<-lm(nitrate~abs, data=stan) #make regression through stand curve data
  Nset$nitrate<-predict(mod, Nset)
  nitratesamples<-rbind(nitratesamples, Nset)
}

nitratesamples<-filter(nitratesamples, !is.na(nitratesamples$nitrate))


