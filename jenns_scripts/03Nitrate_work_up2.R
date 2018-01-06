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

#according to notes, remove spp 1 , temp 24, 31, 38, after May 18
#according to notes, remove spp 2 and 3 , after June 6
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

#make subset of data in which initial concentrations are removed
edata<-subset(cdata, cdata$day!="0")
edata<-subset(edata, edata$species!="1"|edata$day!="7")
dim(edata)

edata$species<-as.factor(edata$species)
edata$temp<-as.factor(edata$temp)

#plot nitrate by temp and species along day and fit equilibrium flat line model
par(mfcol=c(6,4))
par(mar=c(3,3,0,0), oma=c(1,1,0.5, 0.5))
results<-data.frame(species=c(rep(1,6),rep(2,6),rep(3,6),rep(4,6)), temp=c(rep(c(3, 10, 17, 24, 31, 38), 4)), 
                    nstar=rep(99, 24))
for(j in 1:4){
  for (i in c(3, 10, 17, 24, 31, 38)){
    tempdata<-subset(cdata, cdata$temp==i & cdata$species==j)
    with(tempdata, plot(nitrate~day, ylim=c(0, 12), col=j))
    
    tempdata<-subset(edata, edata$temp==i & edata$species==j)
    model<-lme(nitrate~1, data=tempdata, random=~1|day)
    results$nstar[results$species==j & results$temp==i]<-as.numeric(fixed.effects(model))
    results$nstarSE[results$species==j & results$temp==i]<-as.numeric(sqrt(summary(model)$varFix))
    abline(as.numeric(fixed.effects(model)),0, col=j)
}
}

#plot Rstar as a function of temp
par(mfrow=c(1,1))
with(results, plot(nstar~temp, ylim=c(1, 5), type="n"))
for(i in 1:4){
  with(subset(results, results$species==i), points(nstar~temp, col=i))
  with(subset(results, results$species==i), lines(nstar~temp, col=i))
  with(subset(results, results$species==i), segments(temp, nstar-nstarSE, temp, nstar+nstarSE, col=i))
}



#plot by temp and species along day and fit model
par(mfcol=c(6,4))
par(mar=c(3,3,0,0), oma=c(1,1,0.5, 0.5))
results<-data.frame(species=c(rep(1,6),rep(2,6),rep(3,6),rep(4,6)), temp=c(rep(c(3, 10, 17, 24, 31, 38), 4)))
for(j in 1:4){
for (i in c(3, 10, 17, 24, 31, 38)){
  tryCatch({
  tempdata<-subset(cdata, cdata$temp==i & cdata$species==j)
  with(tempdata, plot(n.fixed~day, ylim=c(0, 12), col=j))

  #model
  a_start = 2
  b_start = 0.05
  y_start = 2.5
  curve <- nls(n.fixed~y+a*exp(-b*day), data=tempdata, 
               start=list(a=a_start,b=b_start, y=y_start), control = list(maxiter = 500))
  
  
  new.data <- data.frame(day = seq(min(tempdata$day),max(tempdata$day),len = 100))
  lines(new.data$day,predict(curve,newdata = new.data))
  results$decay_int[results$species==j & results$temp==i]<-summary(curve)$parameters[3,1]
  results$decay_int_SE[results$species==j & results$temp==i]<-summary(curve)$parameters[3,2]
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  tempdata<-subset(edata, edata$temp==i & edata$species==j)
  model<-lme(nitrate~1, data=tempdata, random=~1|day)
  results$flatline_int[results$species==j & results$temp==i]<-as.numeric(fixed.effects(model))
  results$flatline_int_SE[results$species==j & results$temp==i]<-as.numeric(sqrt(summary(model)$varFix))
  abline(as.numeric(fixed.effects(model)),0, col=j)
}
}


  for (i in c(3, 10, 38)){
tempdata<-subset(cdata, cdata$species==1 & cdata$temp==i & cdata$day>5)
with(tempdata, plot(n.fixed~day, ylim=c(0, 12), col=j))

#model
a_start = 2
b_start = 0.05
y_start = 2.5
curve <- nls(n.fixed~y+a*exp(-b*day), data=tempdata, 
             start=list(a=a_start,b=b_start, y=y_start), control = list(maxiter = 500))


new.data <- data.frame(day = seq(min(tempdata$day),max(tempdata$day),len = 100))
lines(new.data$day,predict(curve,newdata = new.data))
results$delay_decay_int[results$species==1 & results$temp==i]<-summary(curve)$parameters[3,1]
results$delay_decay_int_SE[results$species==1 & results$temp==i]<-summary(curve)$parameters[3,2]
}

curve(y_start+(a_start)*exp((-b_start)*x), 0, 50, add=TRUE)

#plot Rstar as a function of temp
par(mfrow=c(1,1))
with(results, plot(flatline_int~temp, type="n", ylim=c(2, 5)))
for(i in 1:4){
  with(subset(results, results$species==i), points(flatline_int~temp, col=i))
  with(subset(results, results$species==i), lines(flatline_int~temp, col=i))
  with(subset(results, results$species==i), segments(
    temp, flatline_int-flatline_int_SE, temp, flatline_int+flatline_int_SE, col=i))
}



#what was initial [N]?
(0.05-coef(st.2806)[1])/coef(st.2806)[2] #4.4 uM. Not what I thought. all medium samples are coming out low like this.

#temporary
#subset data to just the data for which I made a standard curve
#data<-subset(data, data$date.measured<"2017-07-02"| data$date.measured>"2017-07-11")

#temporary
#make points from May 2 and 5 lower by abs = 0.015
#data$meanabscorr<-ifelse(data$date.sampled=="2017-05-02" | data$date.sampled=="2017-05-05", data$meanabs-0.015, data$meanabs)
#make points from May 16 and 19 higher by abs = 0.015
#data$meanabscorr<-ifelse(data$date.sampled=="2017-05-19", data$meanabscorr+0.015, data$meanabscorr)

#
#temporary
#
#data$meanabs<-data$meanabscorr #replace raw data with "corrected data" for offset in some sampling dates. This is 


#temporary
#subset data to remove points above 7.5 -> errors
plot(data$nitrate, col=ifelse(data$nitrate>7.5, 2, 1))
data<-subset(data, data$nitrate<6)
dim(data)

sp1<-subset(data, data$species==1)
sp2<-subset(data, data$species==2)
sp3<-subset(data, data$species==3)
sp4<-subset(data, data$species==4)

sp1$day<-sp1$date-min(sp1$date)
sp2$day<-sp2$date-min(sp2$date)
sp3$day<-sp3$date-min(sp3$date)
sp4$day<-sp4$date-min(sp4$date)

# have a look at initial values
initsp1<-subset(sp1, sp1$day==0)
initsp2<-subset(sp2, sp1$day==0)
initsp3<-subset(sp3, sp1$day==0)
initsp4<-subset(sp4, sp4$day==0)
init<-rbind(initsp1, initsp2, initsp3, initsp4)

initmean<-mean(init$nitrate, na.rm=T)
initsd<-sd(init$nitrate, na.rm=T)

#subset each species' data to remove day=0 data
sp1<-subset(sp1, sp1$day!=0)
sp2<-subset(sp2, sp2$day!=0)
sp3<-subset(sp3, sp3$day!=0)
sp4<-subset(sp4, sp4$day!=0)

datadays<-rbind(sp1, sp2, sp3, sp4)

#####################################
#model
#lme of nitrate~day+temp*species, random=date.measured
###########################################
datadays$machineday<-as.factor(datadays$date.measured)
datadays$species<-as.factor(datadays$species)
datadays$temp<-as.factor(datadays$temp)
mod<-lme(nitrate~day+temp*species, random=~1|machineday, data=datadays)
summary(mod)
coef(mod)[1,2]

par(mfrow=c(2,2))
plot(nitrate~day, data=datadays, type="n", main="species 1")
abline(mean(coef(mod)[,1]),coef(mod)[1,2], col="blue")
abline(mean(coef(mod)[,1])+coef(mod)[1,3], coef(mod)[1,2], col="turquoise")
abline(mean(coef(mod)[,1])+coef(mod)[1,4], coef(mod)[1,2], col="green")
abline(mean(coef(mod)[,1])+coef(mod)[1,5], coef(mod)[1,2], col="yellow")
abline(mean(coef(mod)[,1])+coef(mod)[1,6], coef(mod)[1,2], col="orange")
abline(mean(coef(mod)[,1])+coef(mod)[1,7], coef(mod)[1,2], col="red")

plot(nitrate~day, data=datadays, type="n", main="species 2")
abline(mean(coef(mod)[,1])+coef(mod)[1,8],coef(mod)[1,2], col="blue")
abline(mean(coef(mod)[,1])+coef(mod)[1,3]+coef(mod)[1,8]+coef(mod)[1,11], coef(mod)[1,2], col="turquoise")
abline(mean(coef(mod)[,1])+coef(mod)[1,4]+coef(mod)[1,8]+coef(mod)[1,12], coef(mod)[1,2], col="green")
abline(mean(coef(mod)[,1])+coef(mod)[1,5]+coef(mod)[1,8]+coef(mod)[1,13], coef(mod)[1,2], col="yellow")
abline(mean(coef(mod)[,1])+coef(mod)[1,6]+coef(mod)[1,8]+coef(mod)[1,14], coef(mod)[1,2], col="orange")
abline(mean(coef(mod)[,1])+coef(mod)[1,7]+coef(mod)[1,8]+coef(mod)[1,15], coef(mod)[1,2], col="red")

plot(nitrate~day, data=datadays, type="n", main="species 3")
abline(mean(coef(mod)[,1])+coef(mod)[1,9],coef(mod)[1,2], col="blue")
abline(mean(coef(mod)[,1])+coef(mod)[1,3]+coef(mod)[1,9]+coef(mod)[1,16], coef(mod)[1,2], col="turquoise")
abline(mean(coef(mod)[,1])+coef(mod)[1,4]+coef(mod)[1,9]+coef(mod)[1,17], coef(mod)[1,2], col="green")
abline(mean(coef(mod)[,1])+coef(mod)[1,5]+coef(mod)[1,9]+coef(mod)[1,18], coef(mod)[1,2], col="yellow")
abline(mean(coef(mod)[,1])+coef(mod)[1,6]+coef(mod)[1,9]+coef(mod)[1,19], coef(mod)[1,2], col="orange")
abline(mean(coef(mod)[,1])+coef(mod)[1,7]+coef(mod)[1,9]+coef(mod)[1,20], coef(mod)[1,2], col="red")

plot(nitrate~day, data=datadays, type="n", main="species 4")
abline(mean(coef(mod)[,1])+coef(mod)[1,10],coef(mod)[1,2], col="blue")
abline(mean(coef(mod)[,1])+coef(mod)[1,3]+coef(mod)[1,10]+coef(mod)[1,21], coef(mod)[1,2], col="turquoise")
abline(mean(coef(mod)[,1])+coef(mod)[1,4]+coef(mod)[1,10]+coef(mod)[1,22], coef(mod)[1,2], col="green")
abline(mean(coef(mod)[,1])+coef(mod)[1,5]+coef(mod)[1,10]+coef(mod)[1,23], coef(mod)[1,2], col="yellow")
abline(mean(coef(mod)[,1])+coef(mod)[1,6]+coef(mod)[1,10]+coef(mod)[1,24], coef(mod)[1,2], col="orange")
abline(mean(coef(mod)[,1])+coef(mod)[1,7]+coef(mod)[1,10]+coef(mod)[1,25], coef(mod)[1,2], col="red")

#species1
sp1temp3<-subset(sp1, sp1$temp==3)
sp1summ3 <- plyr::ddply(sp1temp3, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp1summ3 <- rbind(sp1summ3, c(0, initmean, initsd))
sp1summ3$day<-as.numeric(sp1summ3$day)

sp1temp10<-subset(sp1, sp1$temp==10) 
sp1summ10 <- plyr::ddply(sp1temp10, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp1summ10 <- rbind(sp1summ10, c(0, initmean, initsd))
sp1summ10$day<-as.numeric(sp1summ10$day)

sp1temp17<-subset(sp1, sp1$temp==17) 
sp1summ17 <- plyr::ddply(sp1temp17, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp1summ17 <- rbind(sp1summ17, c(0, initmean, initsd))
sp1summ17$day<-as.numeric(sp1summ17$day)

sp1temp24<-subset(sp1, sp1$temp==24) 
sp1summ24 <- plyr::ddply(sp1temp24, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp1summ24 <- rbind(sp1summ24, c(0, initmean, initsd))
sp1summ24$day<-as.numeric(sp1summ24$day)

sp1temp31<-subset(sp1, sp1$temp==31)
sp1summ31 <- plyr::ddply(sp1temp31, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp1summ31 <- rbind(sp1summ31, c(0, initmean, initsd))
sp1summ31$day<-as.numeric(sp1summ31$day)


sp1temp38<-subset(sp1, sp1$temp==38) 
sp1summ38 <- plyr::ddply(sp1temp38, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp1summ38 <- rbind(sp1summ38, c(0, initmean, initsd))
sp1summ38$day<-as.numeric(sp1summ38$day)

#species2
sp2temp3<-subset(sp2, sp2$temp==3)
sp2summ3 <- plyr::ddply(sp2temp3, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp2summ3 <- rbind(sp2summ3, c(0, initmean, initsd))
sp2summ3$day<-as.numeric(sp2summ3$day)

sp2temp10<-subset(sp2, sp2$temp==10) 
sp2summ10 <- plyr::ddply(sp2temp10, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp2summ10 <- rbind(sp2summ10, c(0, initmean, initsd))
sp2summ10$day<-as.numeric(sp2summ10$day)

sp2temp17<-subset(sp2, sp2$temp==17) 
sp2summ17 <- plyr::ddply(sp2temp17, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp2summ17 <- rbind(sp2summ17, c(0, initmean, initsd))
sp2summ17$day<-as.numeric(sp2summ17$day)

sp2temp24<-subset(sp2, sp2$temp==24) 
sp2summ24 <- plyr::ddply(sp2temp24, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp2summ24 <- rbind(sp2summ24, c(0, initmean, initsd))
sp2summ24$day<-as.numeric(sp2summ24$day)

sp2temp31<-subset(sp2, sp2$temp==31)
sp2summ31 <- plyr::ddply(sp2temp31, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp2summ31 <- rbind(sp2summ31, c(0, initmean, initsd))
sp2summ31$day<-as.numeric(sp2summ31$day)


sp2temp38<-subset(sp2, sp2$temp==38) 
sp2summ38 <- plyr::ddply(sp2temp38, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp2summ38 <- rbind(sp2summ38, c(0, initmean, initsd))
sp2summ38$day<-as.numeric(sp2summ38$day)

####species 3
sp3temp3<-subset(sp3, sp3$temp==3)
sp3summ3 <- plyr::ddply(sp3temp3, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp3summ3 <- rbind(sp3summ3, c(0, initmean, initsd))
sp3summ3$day<-as.numeric(sp3summ3$day)

sp3temp10<-subset(sp3, sp3$temp==10) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp3summ10 <- plyr::ddply(sp3temp10, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp3summ10 <- rbind(sp3summ10, c(0, initmean, initsd))
sp3summ10$day<-as.numeric(sp3summ10$day)

sp3temp17<-subset(sp3, sp3$temp==17) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp3summ17 <- plyr::ddply(sp3temp17, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp3summ17 <- rbind(sp3summ17, c(0, initmean, initsd))
sp3summ17$day<-as.numeric(sp3summ17$day)

sp3temp24<-subset(sp3, sp3$temp==24) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp3summ24 <- plyr::ddply(sp3temp24, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp3summ24 <- rbind(sp3summ24, c(0, initmean, initsd))
sp3summ24$day<-as.numeric(sp3summ24$day)

sp3temp31<-subset(sp3, sp3$temp==31) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp3summ31 <- plyr::ddply(sp3temp31, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp3summ31 <- rbind(sp3summ31, c(0, initmean, initsd))
sp3summ31$day<-as.numeric(sp3summ31$day)

sp3temp38<-subset(sp3, sp3$temp==38) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp3summ38 <- plyr::ddply(sp3temp38, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp3summ38 <- rbind(sp3summ38, c(0, initmean, initsd))
sp3summ38$day<-as.numeric(sp3summ38$day)


####species 4
sp4temp3<-subset(sp4, sp4$temp==3)
sp4summ3 <- plyr::ddply(sp4temp3, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp4summ3 <- rbind(sp4summ3, c(0, initmean, initsd))
sp4summ3$day<-as.numeric(sp4summ3$day)

sp4temp10<-subset(sp4, sp4$temp==10) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp4summ10 <- plyr::ddply(sp4temp10, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp4summ10 <- rbind(sp4summ10, c(0, initmean, initsd))
sp4summ10$day<-as.numeric(sp4summ10$day)

sp4temp17<-subset(sp4, sp4$temp==17) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp4summ17 <- plyr::ddply(sp4temp17, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp4summ17 <- rbind(sp4summ17, c(0, initmean, initsd))
sp4summ17$day<-as.numeric(sp4summ17$day)

sp4temp24<-subset(sp4, sp4$temp==24) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp4summ24 <- plyr::ddply(sp4temp24, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp4summ24 <- rbind(sp4summ24, c(0, initmean, initsd))
sp4summ24$day<-as.numeric(sp4summ24$day)

sp4temp31<-subset(sp4, sp4$temp==31) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp4summ31 <- plyr::ddply(sp4temp31, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp4summ31 <- rbind(sp4summ31, c(0, initmean, initsd))
sp4summ31$day<-as.numeric(sp4summ31$day)

sp4temp38<-subset(sp4, sp4$temp==38) #I am here - trying to subset to remove rows with abs.1 and abs.2 are both NA
sp4summ38 <- plyr::ddply(sp4temp38, "day", plyr::summarise, mean = mean(nitrate), sd = sd(nitrate))
sp4summ38 <- rbind(sp4summ38, c(0, initmean, initsd))
sp4summ38$day<-as.numeric(sp4summ38$day)

#
#
#
#
