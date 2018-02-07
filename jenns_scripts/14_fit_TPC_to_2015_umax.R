#goal: fit Schoolfield TPC to umax data from 2015 
#to get an estimate of umax for 200 temps between 3 and 38, 100 runs from 4 species
#currently just TT

library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

#read in data
TT_ks_umax<-cbind(read_csv("data-processed/TT_ks_umax_boot.csv"), Species="TT")
CS_ks_umax<-cbind(read_csv("data-processed/CS_ks_umax_boot.csv"), Species="CS")
AC_ks_umax<-cbind(read_csv("data-processed/AC_ks_umax_boot.csv"), Species="AC")
CH_ks_umax<-cbind(read_csv("data-processed/CH_ks_umax_boot.csv"), Species="CH")
all_ks_umax<-rbind(TT_ks_umax, CS_ks_umax, AC_ks_umax, CH_ks_umax)
all_umax<- all_ks_umax %>%
  filter(term=="umax") 


# get data in order -------------------------------------------------------
#make a subset of first run, all temperatures, change temp to Kelvins
current_dataset <- all_umax %>% 
  rename(umax=estimate) %>%
  filter(run==1, Species=="TT") %>%
  select(umax, Temperature) %>% 
  mutate(K = Temperature + 273.15) %>% 
  rename(OriginalTraitValue = umax) %>% 
  select(-Temperature)
current_dataset$OriginalTraitValue[current_dataset$OriginalTraitValue == 0] <- 1 #change zeros to 1s

## If there are negative values, substract the minimum value
MinVal <- NA
if (min(current_dataset$OriginalTraitValue)<=0){
  MinVal <- min(current_dataset$OriginalTraitValue)
  current_dataset$OriginalTraitValue <-current_dataset$OriginalTraitValue - MinVal
  current_dataset <-current_dataset[-which(current_dataset$OriginalTraitValue==0),]}



#### assign Tref as GlobalEnv
# T_ref is the standardization temperature (in K). 
# This needs to be any value below the peak of the curve.
assign("Tref", 285.15, envir = .GlobalEnv) 



# Estimate STARTING VALUES for the nls ------------------------------------

GetE <- function(tmp, rate, T.p, k=8.62e-5)
{
  # Estimate starting value for E, taking linear regression using the rise part
  # of the curve only.
  # ~~~ Parameters ~~~
  # tmp  : temperature data (in K).
  # rate : rate data corresponding to temperature above.
  # T.p  : temperature at which rate peaks, used as a cutoff point.
  # k    : Boltzmann constant.
  
  tmp.w <- which(tmp <= T.p)
  if (length(tmp.w) > 1)
  {
    m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp[tmp.w]))))
    return(abs(summary(m)$coefficients[2, 1]))
  } else
  {
    return(0.6)
  }
}

GetB0 <- function(tmp, rate)
{
  # Estimate starting value for the normalising constant.
  # ~~~ Parameters ~~~
  # tmp   : temperature data (in K).
  # rate  : rate data corresponding to temperature above.
  # T.ref : estimate normalising constant at this temperature (in K).
  
  if (min(tmp,na.rm=TRUE) > Tref)
  {
    return(log(min(rate[1],na.rm=TRUE)))
  } else
  {
    return(log(max(rate[which(tmp <= Tref)],na.rm=TRUE)))
  }
}


GetTpk <- function(tmp, rate)
{
  # Temperature at which the rate is maximised (estimate of T.peak).
  # ~~~ Parameters ~~~
  # tmp  : Temperature data (in K).
  # rate : Rate data corresponding to temperature above.
  
  return(max(tmp[which.max(rate)]))
}




# Schoolfield fitting -----------------------------------------------------

Schoolfield <- function(B0, E, E_D, T_h, temp) {
  
  # Boltzmann's constant. Units imply that E and E_D are in eV.
  k <- 8.62e-5
  
  # B0 is the normalization constant.    
  # E is the activation energy.
  # E_D is the de-activation energy.    
  # T_h is the temperature at which the rate-limiting enzyme 
  # is 50% active and 50% denatured due to high temperature.
  
  #     return(B0 - E/k * (1/temp - 1/Tref) - log(1 + exp((E_D/k) * (1/T_h - 1/temp)))) #Schoolfied model in original form (not with T_pk as an explicit parameter)
  
  return(B0 + log(exp((-E/k) * ((1/temp) - (1/Tref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_h - 1/temp))))) ## T_pk as an explicit parameter. FITS BETTER
  
}


B0_sch <- c()
E_sch <- c()
E_D_sch <- c()	
T_h_sch <- c()
T_pk_sch <- c()
P_pk_sch <- c()
AIC_sch <- c()
r_sq_sch <- c()

T.h.st  <- GetTpk(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)
E.st    <- GetE(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue, T.p=T.h.st)
B.st <- GetB0(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)


#
#
#
#
#

# fitting schoolfield model to data-------------------------------------------------------------
#start with just TT for now but could do for each species as well
allruns_dataset <- all_umax %>% 
  rename(umax=estimate) %>%
  filter(Species=="TT") %>%
  select(umax, Temperature, run) %>% 
  mutate(K = Temperature + 273.15) %>% 
  rename(OriginalTraitValue = umax)

with(allruns_dataset, plot(OriginalTraitValue~Temperature, 
                           xlim=c(3, 38), ylim=c(0, max(OriginalTraitValue)*1.1)))

for(i in 1:100){
  current_dataset<- allruns_dataset %>%
    group_by(Temperature) %>%
    sample_n(., 1)

schoolfield_nls <- NA
schoolfield_nls <- nlsLM(
  log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
  start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
  lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
  upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=273.15+150),
  data=current_dataset, control=list(minFactor=1 / 2^16, maxiter=1024))

if(!is.na(schoolfield_nls[1])) 
{ 
  
  # Collect the parameter estimates by adding the currect estimate to the vector of previous estimates
  E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
  E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
  T_h_sch <- c(T_h_sch, coef(schoolfield_nls)["T_h"])
  AIC_sch<- c(AIC_sch, AIC(schoolfield_nls))
  
  # Calculate the R squared value as: 1 - (rss/tss)
  rss <- sum((exp(predict(schoolfield_nls)) - 
                current_dataset$OriginalTraitValue)^2, 
             na.rm = TRUE)
  tss <- sum(
    (current_dataset$OriginalTraitValue - 
       mean(current_dataset$OriginalTraitValue, na.rm = TRUE))^2, 
    na.rm = TRUE)
  
  if ( tss != 0 )
  {
    r_sq_sch <- c(r_sq_sch, 1 - (rss/tss))
  } else
  {
    r_sq_sch <- c(r_sq_sch, 1)
  }
  
  # Calculate the peak of the curve and its 
  # corresponding temperature value.
  curr_prediction <- predict(schoolfield_nls)
  for (j in 1:length(curr_prediction))
  {
    # If we found the maximum performance, exit the loop.
    if (curr_prediction[j] == max(curr_prediction))
    {
      break
    }
  }
  
  T_pk_sch <- c(T_pk_sch, current_dataset$K[j])
}	



##############################
# Plotting Schoolfield's fit #
##############################

# Create a name for the output file using:
#	- the original id number
#   - the species name
#   - the model
# 
# 
# Generate predictions from the model fit...
tmp_temps <- seq(3 + 273.15, 
                 38 + 273.15, length = 200)


tmp_model <- exp(Schoolfield(
  coef(schoolfield_nls)["B0"],
  coef(schoolfield_nls)["E"],
  coef(schoolfield_nls)["E_D"],
  coef(schoolfield_nls)["T_h"],
  tmp_temps
))

ModelToPlotS <- data.frame(
  Temperature = tmp_temps - 273.15, 
  TraitValue = tmp_model, run=i
)


# Prepare the data points of the original values.
DataToPlot <- data.frame(
  Temperature = current_dataset$K - 273.15, 
  TraitValue = current_dataset$OriginalTraitValue
)

with(ModelToPlotS, lines(TraitValue~Temperature))
with(DataToPlot, points(TraitValue~Temperature))

name<-data.frame(umax=ModelToPlotS$TraitValue, Temperature=ModelToPlotS$Temperature, run=i, Species="TT")
assign(paste("Schoolfield_TPC", i, "TT", sep ="_"), name)
}

# bind all of the objects that start with Schoolfield_TPC
TT_Schoolfield_TPC_all<-do.call("rbind", mget(ls(pattern="Schoolfield_TPC")))
dim(TT_Schoolfield_TPC_all)
head(TT_Schoolfield_TPC_all)
write_csv(TT_Schoolfield_TPC_all, "data-processed/TT_Schoolfield_TPC_all.csv")

plot(TT_Schoolfield_TPC_all$umax~TT_Schoolfield_TPC_all$Temperature, type="l")

