#Goal: create functions to fit logistic growth curve to nitrate data

#read in data
cdata<-read_csv("data-processed/processed_nitrate_data.csv")
names(cdata)

#read in libraries
library(simecol)
library(tidyverse)

#figure out the mean starting condition
N_init<-4.41639

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
  with(observeddata, points(N~obstime, ylim=c(0, 12)))
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



par(mfcol=c(4,6))

#for species 1, temp 1, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.03, K = 1)
UpperBound <- c(r = 2, K = 5) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,75), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==3 & cdata$species==1 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results031<-data.frame(r=resultsr*(-1), K=resultsk, temp=3, species=1)

#for species 2, temp 1, use:
Parameters <- c(r = 0.5, K = 0.01)
LowerBound <- c(r = 0.06, K = 2)
UpperBound <- c(r = 2, K = 5) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==3 & cdata$species==2 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}
results032<-data.frame(r=resultsr*c(-1, 1, -1, 1, 1), K=resultsk, temp=3, species=2)

#for species 3, temp 1, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.6, K = 2)
UpperBound <- c(r = 2, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==3 & cdata$species==3 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results033<-data.frame(r=resultsr*c(-1, 1, -1, 1, 1), K=resultsk, temp=3, species=3)

#for species 4, temp 1, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 2)
UpperBound <- c(r = 2, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==3 & cdata$species==4 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results034<-data.frame(r=resultsr*(-1), K=resultsk, temp=3, species=4)

#for species 1, temp=10, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==10 & cdata$species==1 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results101<-data.frame(r=resultsr*(-1), K=resultsk, temp=10, species=1)

#for species 2, temp=10, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==10 & cdata$species==2 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results102<-data.frame(r=resultsr*(-1), K=resultsk, temp=10, species=2)
#what is going on with rep 4?

#for species 3, temp=10, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==10 & cdata$species==3 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results103<-data.frame(r=resultsr*(-1), K=resultsk, temp=10, species=3)

#for species 4, temp=10, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==10 & cdata$species==4 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results104<-data.frame(r=resultsr*(-1), K=resultsk, temp=10, species=4)

#for species 1, temp=17, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==17 & cdata$species==1 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results171<-data.frame(r=resultsr*(-1), K=resultsk, temp=17, species=1)

#for species 2, temp=17, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.09, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==17 & cdata$species==2 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results172<-data.frame(r=resultsr*(-1), K=resultsk, temp=17, species=2)

#for species 3, temp=17, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==17 & cdata$species==3 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results173<-data.frame(r=resultsr*(-1), K=resultsk, temp=17, species=3)

#for species 4, temp=17, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==17 & cdata$species==4 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results174<-data.frame(r=resultsr*(-1), K=resultsk, temp=17, species=4)

#for species 1, temp=24, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==24 & cdata$species==1 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results241<-data.frame(r=resultsr*(-1), K=resultsk, temp=24, species=1)

#for species 2, temp=24, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==24 & cdata$species==2 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results242<-data.frame(r=resultsr*(-1), K=resultsk, temp=24, species=2)
#consider dropping rep 5 because it is influenced by very high N values unrealistic within the experiment

#for species 3, temp=24, use:
Parameters <- c(r = 3, K = 5)
LowerBound <- c(r = 1, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.01 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==24 & cdata$species==3 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results243<-data.frame(r=resultsr*(-1), K=resultsk, temp=24, species=3)

#for species 4, temp=24, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==24 & cdata$species==4 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results244<-data.frame(r=resultsr*(-1), K=resultsk, temp=24, species=4)

#for species 1, temp=31, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==31 & cdata$species==1 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results311<-data.frame(r=resultsr*(-1), K=resultsk, temp=31, species=1)

#for species 2, temp=31, use:
Parameters <- c(r = 0.5, K = 5)
LowerBound <- c(r = 0.06, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==31 & cdata$species==2 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results312<-data.frame(r=resultsr*(-1), K=resultsk, temp=31, species=2)

#for species 3, temp=31, use:
Parameters <- c(r = 2, K = 5)
LowerBound <- c(r = 0.07, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==31 & cdata$species==3 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results313<-data.frame(r=resultsr*(-1), K=resultsk, temp=31, species=3)

#for species 4, temp=31, use:
Parameters <- c(r = 2, K = 5)
LowerBound <- c(r = 0.07, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==31 & cdata$species==4 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results314<-data.frame(r=resultsr*(-1), K=resultsk, temp=31, species=4)

#for species 1, temp=38, use:
Parameters <- c(r = 2, K = 5)
LowerBound <- c(r = 1, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==38 & cdata$species==1 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results381<-data.frame(r=resultsr*(-1), K=resultsk, temp=38, species=1)
#fixed r at 1 because otherwise the r was sometimes very low, suggesting N updake occured
#slowly in the system, but we know that cells dies early on so most of action should be at beginning.

#for species 2, temp=38, use:
Parameters <- c(r = 2, K = 5)
LowerBound <- c(r = 1, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==38 & cdata$species==2 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results382<-data.frame(r=resultsr*(-1), K=resultsk, temp=38, species=2)
#fixed r at 1 because otherwise the r was sometimes very low, suggesting N updake occured
#slowly in the system, but we know that cells dies early on so most of action should be at beginning.

#for species 3, temp=38, use:
Parameters <- c(r = 2, K = 5)
LowerBound <- c(r = 1, K = 1)
UpperBound <- c(r = 4, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==38 & cdata$species==3 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results383<-data.frame(r=resultsr*(-1), K=resultsk, temp=38, species=3)
#fixed r at 1 because otherwise the r was sometimes very low, suggesting N updake occured
#slowly in the system, but we know that cells dies early on so most of action should be at beginning.

#for species 4, temp=38, use:
Parameters <- c(r = 2, K = 5)
LowerBound <- c(r = 1, K = 1)
UpperBound <- c(r = 5, K = 6) 
ParamScaling <- 0.001 / UpperBound

plot(c(1,50), c(0,12), type="n")
for(j in 1:length(unique(cdata$rep))){
  curvedata<-subset(cdata, cdata$temp==38 & cdata$species==4 & cdata$rep==j)
  fitgrowth<-plotsinglefit(curvedata)
  resultsr[j]<-fitgrowth$r[1]
  resultsk[j]<-fitgrowth$K[1]
}

results384<-data.frame(r=resultsr*c(-1, 1, 1, 1, -1), K=resultsk, temp=38, species=4)
#fixed r at 1 because otherwise the r was sometimes very low, suggesting N updake occured
#slowly in the system, but we know that cells dies early on so most of action should be at beginning.

#bring results together

results03<-rbind(results031, results032, results033, results034) %>%
  mutate(rep=rep(1:5, 4))
results10<-rbind(results101, results102, results103, results104) %>%
  mutate(rep=rep(1:5, 4))
results17<-rbind(results171, results172, results173, results174) %>%
  mutate(rep=rep(1:5, 4))
results24<-rbind(results241, results242, results243, results244) %>%
  mutate(rep=rep(1:5, 4))
results31<-rbind(results311, results312, results313, results314) %>%
  mutate(rep=rep(1:5, 4))
results38<-rbind(results381, results382, results383, results384) %>%
  mutate(rep=rep(1:5, 4))

nitrate_decay_results<-rbind(results03, results10, results17, results24, results31, results38)

#add species names
nitrate_decay_results_with_species<-nitrate_decay_results %>%
  mutate(Species=ifelse(species==1, "TT", ifelse(species==2, "CS", ifelse(species==3, "AC", ifelse(species==4, "CH", NA)))))
write.csv(nitrate_decay_results_with_species, file="data-processed/nitrate_decay_rep_r-star.csv")

#
#
#
nitrate_decay_results=read.csv(file="data-processed/nitrate_decay_rep_r-star.csv")

nitrate_decay_results %>%
  filter(!temp==38 | !Species=="TT")%>%
  filter(!temp==38 | !Species=="CS")%>%
  filter(!temp==38 | !Species=="AC")%>%
  ggplot(aes(x = temp, y = K, color=as.factor(Species))) + geom_point(size = 2) +
  theme_bw() + facet_grid(~Species) + geom_smooth(method = 'loess')
ggsave("rstar_with_temp_2017.pdf", width=7, height=2)


par(mfrow=c(1,1))
with(nitrate_decay_results, plot(K~temp, type="n", las=1))
for(i in 1:4){
  with(subset(nitrate_decay_results, nitrate_decay_results$species==i), points(K~temp, col=i))
  with(subset(nitrate_decay_results, nitrate_decay_results$species==i), lines(K~temp, col=i))  
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
