#Goal: create functions to fit logistic growth

#read in data
Rtemp_all<-read.csv("data-processed/Rtemp_all.csv")


#read in libraries
library(simecol)
library(tidyverse)


#fix errors in temperature treatment names
Rtemp_all$temperature[Rtemp_all$temperature=="1"]<-"3"
Rtemp_all$temperature[Rtemp_all$temperature=="30"]<-"31"
Rtemp_all$P<- Rtemp_all$cell_density


#troubleshoot some weird looking data
subset(Rtemp_all, Rtemp_all$temperature=="31" & Rtemp_all$species==1 & Rtemp_all$cell_density>2000)

#figure out the mean starting condition and put that into model in previous script
mean(subset(Rtemp_all, Rtemp_all$time_since_innoc_days<1)$P)

# plot just the data ------------------
par(mfrow=c(2,3))
for(i in 1:6){
  curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[i] & Rtemp_all$species==2)
  with(curvedata, plot((cell_density)~time_since_innoc_days, col=as.numeric(temperature), 
                       main=unique(Rtemp_all$temperature)[i], ylim=c(0, 30000)))
}
names(Rtemp_all)

# set up models -----------------------------------------------------------

Parameters <- c(r = 0.05, K = 25000)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.015, K = 0)
UpperBound <- c(r = 2, K = 50000) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

CRmodel <- new("odeModel",
               main = function (time, init, parms) {
                 with(as.list(c(init, parms)), {
                   dp <-  r * P * (1 - (P / K))
                   list(c(dp))
                 })
               },
               parms = Parameters, # Trying Joey's empirically estimated parameters here.
               times = c(from = 0, to = 120, by = 1), # the time interval over which the model will be simulated.
               init = c(P = 1489),
               solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
)

fittedparms <- c("r", "K") # for assigning fitted parameter values to fittedCRmodel

# fitting function -----------------------------------------------------------
controlfit <- function(curvedata){
  init(CRmodel) <- c(P = 1489) # Set initial model conditions to the biovolume taken from the first measurement day
  obstime <- curvedata$time_since_innoc_days # The X values of the observed data points we are fitting our model to
  yobs <- select(curvedata, P) # The Y values of the observed data points we are fitting our model to
  
  
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
  init(CRmodel) <- c(P = 1489) # Set initial model conditions to the biovolume taken from the first measurement day
  obstime <- curvedata$time_since_innoc_days # The X values of the observed data points we are fitting our model to
  yobs <- select(curvedata, P) # The Y values of the observed data points we are fitting our model to
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
  
  r <- coef(fittedCRmodel)[1]
  K <- coef(fittedCRmodel)[2]
  output <- data.frame(r, K)
  return(output)
  
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
  with(observeddata, plot(P~obstime, ylim=c(0, max(P)*2)))
  with(simulateddata, lines(P~time, col=2))
}

# fit and plot ------------------------------------------------------------
Parameters <- c(r = 0.5, K = 25000)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0, K = 10 ^ 2)
UpperBound <- c(r = 1, K = 40000) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 0.001 / UpperBound

#play with subsets of the data
Rtemp_all$P<- Rtemp_all$cell_density
Rtemp_all$temperature<-as.numeric(Rtemp_all$temperature)
curvedata<-subset(Rtemp_all, Rtemp_all$temperature==unique(Rtemp_all$temperature)[1] & Rtemp_all$species==4)

controlfit(curvedata)
plotsinglefit(curvedata)


