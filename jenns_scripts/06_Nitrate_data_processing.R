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

#convert date into day of experiment
sp1<-subset(data, data$species==1)
sp2<-subset(data, data$species==2)
sp3<-subset(data, data$species==3)
sp4<-subset(data, data$species==4)

sp1$day<-sp1$date-min(sp1$date)
sp2$day<-sp2$date-min(sp2$date)
sp3$day<-sp3$date-min(sp3$date)
sp4$day<-sp4$date-min(sp4$date)

data<-rbind(sp1, sp2, sp3, sp4)
data$day<-as.numeric(data$day)


#remove visual outliers (nitrate values abherrently high, contaminated samples, low volume N assays that were abherrent)
#but leave all day 0
levels(data$error_code)
day0<-subset(data, day==0)
aboveday0<-subset(data, day>0)
caboveday0<-subset(aboveday0, !aboveday0$error_code %in% c("contaminated", "vis_outlier", "remove_low_vol", "discard"))
cdata<-rbind(day0, caboveday0)

#how many vis outliers?
dim(subset(aboveday0, aboveday0$error_code %in% 
             c("vis_outlier")))/dim(data)
test<-subset(aboveday0, aboveday0$error_code %in% 
               c("vis_outlier"))
par(mfrow=c(1,1))
plot(data$absorbance.1, col=ifelse(data$error_code=="vis_outlier", 2, 1))


#take average of 2 absorbances
cdata$meanabs<-(cdata$absorbance.1+cdata$absorbance.2)/2
#remove NAs
#cdata<-subset(cdata, complete.cases(cdata$species))

#convert absorbance to nitrate
#for data before July 2nd, use standard curve from 280617; 
#for data from after July 2nd, use standard curve from 110717 
#(or, consider interpolating standard curve for every day between standard curve creation)


#cdata$nitrate<- all_st_curves_model$estimate[1] + cdata$meanabs*all_st_curves_model$estimate[2]


#use the predict function to estimate mean and standard error of nitrate from standard curve models
#the models:
mod_17.06.26
mod_17.06.27
mod_17.06.28
mod_17.07.11
mod_17.07.12
mod_17.07.14

par(mfrow=c(1,1))
hist(cdata$date.measured, 20)

#data on 2017-06-26
subset<-cdata %>%
  filter(date.measured<="2017-06-26")
newdataframe<-data.frame(subset$meanabs) %>%
  rename(absorbance=subset.meanabs)
predicted_N_vec<-predict(mod_17.06.26, newdata=newdataframe)
predicted_N_SE_vec<-predict(mod_17.06.26, newdata=newdataframe, se.fit=TRUE)$se.fit
subset26<-subset %>%
  mutate(predicted_N=predicted_N_vec, predicted_N_SE=predicted_N_SE_vec)


#data on 2017-06-27
subset<-cdata %>%
  filter(date.measured=="2017-06-27")
newdataframe<-data.frame(subset$meanabs) %>%
  rename(absorbance=subset.meanabs)
predicted_N_vec<-predict(mod_17.06.27, newdata=newdataframe)
predicted_N_SE_vec<-predict(mod_17.06.27, newdata=newdataframe, se.fit=TRUE)$se.fit
subset27<-subset %>%
  mutate(predicted_N=predicted_N_vec, predicted_N_SE=predicted_N_SE_vec)
subset27$date.measured

#data on 2017-06-28 up to 2017-07-02
subset<-cdata %>%
  filter(date.measured>="2017-06-27" & date.measured>="2017-07-02")
newdataframe<-data.frame(subset$meanabs) %>%
  rename(absorbance=subset.meanabs)
predicted_N_vec<-predict(mod_17.06.28, newdata=newdataframe)
predicted_N_SE_vec<-predict(mod_17.06.28, newdata=newdataframe, se.fit=TRUE)$se.fit
subset28to02<-subset %>%
  mutate(predicted_N=predicted_N_vec, predicted_N_SE=predicted_N_SE_vec)

#data on 2017-06-28 up to and including 2017-07-02 (but no data past 17 06 30)
subset<-cdata %>%
  filter(date.measured>="2017-06-27" & date.measured<="2017-07-02")
newdataframe<-data.frame(subset$meanabs) %>%
  rename(absorbance=subset.meanabs)
predicted_N_vec<-predict(mod_17.06.28, newdata=newdataframe)
predicted_N_SE_vec<-predict(mod_17.06.28, newdata=newdataframe, se.fit=TRUE)$se.fit
subset28to02<-subset %>%
  mutate(predicted_N=predicted_N_vec, predicted_N_SE=predicted_N_SE_vec)


#data on 2017-07-03 up to and including 2017-07-11 
subset<-cdata %>%
  filter(date.measured>"2017-07-03" & date.measured<="2017-07-11")
newdataframe<-data.frame(subset$meanabs) %>%
  rename(absorbance=subset.meanabs)
predicted_N_vec<-predict(mod_17.07.11, newdata=newdataframe)
predicted_N_SE_vec<-predict(mod_17.07.11, newdata=newdataframe, se.fit=TRUE)$se.fit
subset03to11<-subset %>%
  mutate(predicted_N=predicted_N_vec, predicted_N_SE=predicted_N_SE_vec)

#data on 2017-07-12
subset<-cdata %>%
  filter(date.measured=="2017-07-12")
newdataframe<-data.frame(subset$meanabs) %>%
  rename(absorbance=subset.meanabs)
predicted_N_vec<-predict(mod_17.07.12, newdata=newdataframe)
predicted_N_SE_vec<-predict(mod_17.07.12, newdata=newdataframe, se.fit=TRUE)$se.fit
subset12<-subset %>%
  mutate(predicted_N=predicted_N_vec, predicted_N_SE=predicted_N_SE_vec)

#data on 2017-07-13 ad higher
subset<-cdata %>%
  filter(date.measured>="2017-07-13")
newdataframe<-data.frame(subset$meanabs) %>%
  rename(absorbance=subset.meanabs)
predicted_N_vec<-predict(mod_17.07.12, newdata=newdataframe)
predicted_N_SE_vec<-predict(mod_17.07.12, newdata=newdataframe, se.fit=TRUE)$se.fit
subset13up<-subset %>%
  mutate(predicted_N=predicted_N_vec, predicted_N_SE=predicted_N_SE_vec)

#bring these all together
cdata_nitrate<-rbind(subset26, subset27, subset28to02, subset03to11, subset12, subset13up)

#cdata$nitrate<-ifelse(cdata$date.measured<"2017-07-02", 
#                     (cdata$meanabs-coef(st.2806)[1])/coef(st.2806)[2], 
#                     ifelse(cdata$date.measured<"2017-07-13" & cdata$date.measured>"2017-07-02", 
#                            (cdata$meanabs-coef(st.1207)[1])/coef(st.1207)[2], 
#                            (cdata$meanabs-coef(st.1407)[1])/coef(st.1407)[2]))

plot(cdata_nitrate$meanabs~cdata_nitrate$predicted_N)


cdata<-cdata_nitrate
#plot data
par(mfrow=c(1,2))
with (cdata, plot(meanabs~date))
with (cdata, plot(predicted_N~date))

#
#
#
#

#according to notes, remove highest temp of species 1 before May 15
sp1<-subset(cdata, cdata$species==1)
sp2<-subset(cdata, cdata$species==2)
sp3<-subset(cdata, cdata$species==3)
sp4<-subset(cdata, cdata$species==4)
sp1c<-subset(sp1, sp1$temp<35 | sp1$date<"2017-05-15")
cdata<-rbind(sp1c, sp2, sp3, sp4)

#set the initial concentrations
day0<-subset(cdata, cdata$day=="0")
day0$n.fixed<-mean(subset(day0$predicted_N, complete.cases(day0$predicted_N)))
daynot0<-subset(cdata, cdata$day!="0")
daynot0$n.fixed<-daynot0$predicted_N
cdata<-rbind(day0, daynot0)

write_csv(cdata, "data-processed/processed_nitrate_data_SE.csv")
