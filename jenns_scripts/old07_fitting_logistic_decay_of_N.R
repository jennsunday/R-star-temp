#Goal: create functions to fit logistic growth curve to nitrate data

#get cdata from Nitrate work up file
names(cdata)

#read in libraries
library(simecol)
library(tidyverse)

#figure out the mean starting condition
N_init<-mean(day0$nitrate)


#define N
#note: "n.fixed" is nitrate column when day 0 N is all fixed; "nitrate" is nitrate column when day 0 are raw values
cdata$N<-cdata$n.fixed

# plot just the data ------------------ 
par(mfcol=c(6,4))
par(mar=c(3,3,0,0), oma=c(1,1,0.5, 0.5))
for(j in 1:4){
  for (i in c(3, 10, 17, 24, 31, 38)){
    tempdata<-subset(cdata, cdata$temp==i & cdata$species==j)
    with(tempdata, plot(N~day, ylim=c(0, 12), col=j)) 
  }
}




# set up models -----------------------------------------------------------

Parameters <- c(r = 0.05, K = 0.01)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.001, K = 0.01)
UpperBound <- c(r = 2, K = 5) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

CRmodel <- new("odeModel",
               main = function (time, init, parms) {
                 with(as.list(c(init, parms)), {
                   dp <-  r * N * (1 - (N / K))
                   list(c(dp))
                 })
               },
               parms = Parameters, # Trying Joey's empirically estimated parameters here.
               times = c(from = 0, to = 72, by = 1), # the time interval over which the model will be simulated.
               init = c(N = N_init),
               solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
)

fittedparms <- c("r", "K") # for assigning fitted parameter values to fittedCRmodel

# fitting function -----------------------------------------------------------
controlfit <- function(curvedata){
  init(CRmodel) <- c(N = N_init) # Set initial model conditions to the biovolume taken from the first measurement day
  obstime <- curvedata$day # The X values of the observed data points we are fitting our model to
  yobs <- select(curvedata, N) # The Y values of the observed data points we are fitting our model to
  
  
  fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
                               debuglevel = 0, fn = ssqOdeModel,
                               method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
                               control = list(trace = T)
  )
  
  r <- coef(fittedCRmodel)[1]
  K <- coef(fittedCRmodel)[2]
  #ID <- data$ID[1]
  output <- data.frame(r, K)
  return(output)
}



# plot function -----------------------------------------------------------
plotsinglefit <- function(curvedata){
  init(CRmodel) <- c(N = N_init) # Set initial model conditions to the Nitrate taken from the first measurement day
  obstime <- curvedata$day # The X values of the observed data points we are fitting our model to
  yobs <- select(curvedata, N) # The Y values of the observed data points we are fitting our model to
  # parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.
  
  # Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
  # squared differences between the experimental data and our modelled data. This
  # is fairly standard, although alternatives do exist.
  
  # The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
  # "lower" is a vector containing the lower bound constraints
  # for the parameter values. This may need tweaking.
  
  fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
                               debuglevel = 0, fn = ssqOdeModel,
                               method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
                               control = list(trace = T)
  )
  
  # To display the fitted results we need to create a new OdeModel object. Here
  # we duplicate CRmodel and then alter it to use our new fitted parameters.
  plotfittedCRmodel <- CRmodel
  parms(plotfittedCRmodel)[fittedparms] <- coef(fittedCRmodel)
  
  # set model parameters to fitted values and simulate again
  times(plotfittedCRmodel) <- c(from=0, to=120, by=1)
  ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))
  
  # Form observed data into a dataframe; the simulated data are already in a dataframe
  observeddata <- data.frame(obstime, yobs)
  simulateddata <- ysim
  
  # Plot the results of our model fitting.
  with(observeddata, plot(N~obstime, ylim=c(0, 12)))
  with(simulateddata, lines(N~time, col=2))
  
  #output fitted parameters
  r <- coef(fittedCRmodel)[1]
  K <- coef(fittedCRmodel)[2]
  output <- data.frame(r, K)
  return(output)
}


#test----

curvedata<-subset(cdata, cdata$temp==3 & cdata$species==1)
controlfit(curvedata)
plotsinglefit(curvedata)

#plot---
par(mfcol=c(4,6))
resultsr<-1:4
resultsk<-1:4
flatline_int<-1:4
flatline_int_SE<-1:4

#for species 1-4, temp 1, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.05, K = 0.01)
UpperBound <- c(r = 2, K = 5) 
ParamScaling <- 0.001 / UpperBound

for(i in 1:4){
  curvedata<-subset(cdata, cdata$temp==3 & cdata$species==i)
 fitgrowth<-plotsinglefit(curvedata)
  resultsr[i]<-fitgrowth$r[1]
  resultsk[i]<-fitgrowth$K[1]
  
  tempdata<-subset(edata, edata$temp==3 & edata$species==i)
  model<-lme(nitrate~1, data=tempdata, random=~1|day)
  flatline_int[i]<-as.numeric(fixed.effects(model))
  flatline_int_SE[i]<-as.numeric(sqrt(summary(model)$varFix))
  abline(as.numeric(fixed.effects(model)),0, col=3)
}
results3<-data.frame(r=resultsr*(c(1, -1, -1, 1)), K=resultsk, 
                     K_flat=flatline_int, K_flat_SE=flatline_int_SE, temp=rep(3, 4), species=1:4)

#for species 1-4, temp 2, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.05, K = 0.01)
UpperBound <- c(r = 4, K = 5) 
ParamScaling <- 0.001 / UpperBound

for(i in 1:4){
  curvedata<-subset(cdata, cdata$temp==10 & cdata$species==i)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[i]<-fitgrowth$r[1]
  resultsk[i]<-fitgrowth$K[1]
  
  tempdata<-subset(edata, edata$temp==10 & edata$species==i)
  model<-lme(nitrate~1, data=tempdata, random=~1|day)
  flatline_int[i]<-as.numeric(fixed.effects(model))
  flatline_int_SE[i]<-as.numeric(sqrt(summary(model)$varFix))
  abline(as.numeric(fixed.effects(model)),0, col=3)
}
results10<-data.frame(r=resultsr*(c(1, -1, -1, 1)), K=resultsk, 
                      K_flat=flatline_int, K_flat_SE=flatline_int_SE, temp=rep(10, 4), species=1:4)

#for species 1-4, temp 3, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.001, K = 0.01)
UpperBound <- c(r = 4, K = 5) 
ParamScaling <- 0.001 / UpperBound

for(i in 1:4){
  curvedata<-subset(cdata, cdata$temp==17 & cdata$species==i)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[i]<-fitgrowth$r[1]
  resultsk[i]<-fitgrowth$K[1]
  
  tempdata<-subset(edata, edata$temp==17 & edata$species==i)
  model<-lme(nitrate~1, data=tempdata, random=~1|day)
  flatline_int[i]<-as.numeric(fixed.effects(model))
  flatline_int_SE[i]<-as.numeric(sqrt(summary(model)$varFix))
  abline(as.numeric(fixed.effects(model)),0, col=3)
}
results17<-data.frame(r=resultsr*(c(1, -1, -1, 1)), K=resultsk, 
                      K_flat=flatline_int, K_flat_SE=flatline_int_SE, temp=rep(17, 4), species=1:4)

#for species 1-4, temp 4, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.001, K = 0.01)
UpperBound <- c(r = 4, K = 5) 
ParamScaling <- 0.001 / UpperBound

for(i in 1:4){
  curvedata<-subset(cdata, cdata$temp==24 & cdata$species==i)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[i]<-fitgrowth$r[1]
  resultsk[i]<-fitgrowth$K[1]
  
  tempdata<-subset(edata, edata$temp==24 & edata$species==i)
  model<-lme(nitrate~1, data=tempdata, random=~1|day)
  flatline_int[i]<-as.numeric(fixed.effects(model))
  flatline_int_SE[i]<-as.numeric(sqrt(summary(model)$varFix))
  abline(as.numeric(fixed.effects(model)),0, col=3)
}
results24<-data.frame(r=resultsr*(c(1, -1, -1, 1)), K=resultsk, 
                      K_flat=flatline_int, K_flat_SE=flatline_int_SE, temp=rep(24, 4), species=1:4)

#for species 1-4, temp 5, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.5, K = 0.01)
UpperBound <- c(r = 4, K = 5) 
ParamScaling <- 0.001 / UpperBound

for(i in 1:4){
  curvedata<-subset(cdata, cdata$temp==31 & cdata$species==i)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[i]<-fitgrowth$r[1]
  resultsk[i]<-fitgrowth$K[1]
  
  tempdata<-subset(edata, edata$temp==31 & edata$species==i)
  model<-lme(nitrate~1, data=tempdata, random=~1|day)
  flatline_int[i]<-as.numeric(fixed.effects(model))
  flatline_int_SE[i]<-as.numeric(sqrt(summary(model)$varFix))
  abline(as.numeric(fixed.effects(model)),0, col=3)
}
results31<-data.frame(r=resultsr*(c(1, -1, -1, 1)), K=resultsk, 
                      K_flat=flatline_int, K_flat_SE=flatline_int_SE, temp=rep(31, 4), species=1:4)

#for species 1-4, temp 6, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.5, K = 0.01)
UpperBound <- c(r = 4, K = 5) 
ParamScaling <- 0.001 / UpperBound

for(i in 1:4){
  curvedata<-subset(cdata, cdata$temp==38 & cdata$species==i)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[i]<-fitgrowth$r[1]
  resultsk[i]<-fitgrowth$K[1]
  
  tempdata<-subset(edata, edata$temp==38 & edata$species==i)
  model<-lme(nitrate~1, data=tempdata, random=~1|day)
  flatline_int[i]<-as.numeric(fixed.effects(model))
  flatline_int_SE[i]<-as.numeric(sqrt(summary(model)$varFix))
  abline(as.numeric(fixed.effects(model)),0, col=3)
}
results38<-data.frame(r=resultsr*(c(1, -1, -1, 1)), K=resultsk, 
                      K_flat=flatline_int, K_flat_SE=flatline_int_SE, temp=rep(38, 4), species=1:4)
#
#
#
#


results<-rbind(results3, results10, results17, results24, results31, results38)
write.csv(results, file="data-processed/logistic_N_decay_fits_r-star.csv")

par(mfrow=c(1,1))
with(results, plot(K~temp, type="n", las=1))
for(i in 1:4){
  with(subset(results, results$species==i), points(K~temp, col=i))
  with(subset(results, results$species==i), lines(K~temp, col=i))  
  with(subset(results, results$species==i), points(K_flat~temp, col=i, pch=2))
  with(subset(results, results$species==i), segments(temp, K_flat+K_flat_SE, temp, K_flat-K_flat_SE, col=i))
}
mtext("R-star", 2, 3)
mtext("Temperature (°C)", 1, 3)
legend("top", c("TT", "CS", "AC", "Chlamy"), lty=1, col=1:4, pch=1, bty="n")

par(mfrow=c(1,1))
with(results, plot(K~temp, type="n", las=1))
for(i in 1:4){
  with(subset(results, results$species==i), points(K~temp, col=i))
  with(subset(results, results$species==i), lines(K~temp, col=i))  
}
abline(0,0)
mtext("K", 2, 3, las=1)
mtext("Temperature (°C)", 1, 3)
legend("topright", c("TT", "CS", "AC", "Chlamy"), lty=1, col=1:4, pch=1, bty="n")
