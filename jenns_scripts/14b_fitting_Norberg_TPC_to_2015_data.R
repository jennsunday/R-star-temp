### fit TPC to 2015 umax data 
#goal:fit Norberg model to umax data for 1000 bootstraps for each of 4 species
#worried this will take a long time... 
#start with median/mean data so that I am only fitting one curve to each species,
#then reduce the "searching" scope so that the time to fit each model is shorter.


library(bbmle)
library(tidyverse)

## load data
View(TT_ks_umax)
#read in data
TT_ks_umax<-cbind(read_csv("data-processed/TT_ks_umax_boot.csv"), Species="TT")
CS_ks_umax<-cbind(read_csv("data-processed/CS_ks_umax_boot.csv"), Species="CS")
AC_ks_umax<-cbind(read_csv("data-processed/AC_ks_umax_boot.csv"), Species="AC")
CH_ks_umax<-cbind(read_csv("data-processed/CH_ks_umax_boot.csv"), Species="CH")
all_ks_umax<-rbind(TT_ks_umax, CS_ks_umax, AC_ks_umax, CH_ks_umax)
all_umax<- all_ks_umax %>%
  filter(term=="umax") 

#get data ready

dat.full <- all_umax %>% 
  mutate(eachcurve=paste(Species,run)) %>%
  rename(curve.id = eachcurve,
         temperature = Temperature,
         growth.rate = estimate)

### Jenn, skip to line 139 to skip over the running of the fitting (this takes a little while),
### and read in the estimated parameters

## fit TPCs, code and approach adapted from Thomas et al. 2016
nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2)
  res
}

# Create new vector of unique curve.id values
curve.id.list<-unique(dat.full$curve.id)	# careful - this may change the order of the unique id - needs to be fixed

# Create empty vectors to populate with parameter values and trait estimates
z.list<-rep(NA, length(curve.id.list))				#Parameter 'z'
w.list<-rep(NA, length(curve.id.list))				#Parameter 'w', which is the temperature niche width
a.list<-rep(NA, length(curve.id.list))				#Parameter 'a'
b.list<-rep(NA, length(curve.id.list))				#Parameter 'b'
topt.list<-rep(NA, length(curve.id.list))			#Topt, Optimum temperature for growth
maxgrowth.list<-rep(NA, length(curve.id.list))		#Maximum growth rate (i.e. growth rate at Topt)
rsqr.list<-rep(NA, length(curve.id.list))			#R^2 values for all fits
s.list<-rep(NA, length(curve.id.list))				#Error
n.list<-rep(NA, length(curve.id.list))				#Number of growth rate measurements used in the curve

# Loop through all curve.id.list values to estimate parameters for all curves

for(i in 1:length(curve.id.list)){
  print(i)
  
  # Take a subset of the data corressponding to the ith curve.id.list value
  dat<-subset(dat.full,dat.full$curve.id==curve.id.list[i])
  
  # guess starting values for parameters 'z' and 'w'
  z.guess<-mean(dat$temperature[dat$growth.rate==max(dat$growth.rate)])		#starting estimates for 'z'
  w.guess<-diff(range(dat$temperature))										#starting estimates for niche width
  

  #reduced selection based in initial exploration:
  avals<-seq(0,0.2,0.1) 	
  bvals<-seq(0,0.2,0.1) 	

  mod.list<-list()
  AIC.list<-c()
  
  for(ia in 1:length(avals)){
    for(ib in 1:length(bvals)){
      a.guess<-avals[ia]
      b.guess<-bvals[ib]
      res2<-try(fit<-mle2(dat$growth.rate~dnorm(mean=nbcurve(dat$temperature,z=z,w=w,a=a,b=b),sd=s),start=list(z=z.guess,w=w.guess,a=a.guess,b=b.guess,s=0.3),
                          skip.hessian=TRUE,data=dat))
      if(class(res2)!="try-error"){
        mod.list<-append(mod.list,fit)
        AIC.list<-append(AIC.list,AIC(fit))
      }
    }
  }
  
  # Identify the best model from the list and save coefficients and R^2 values
  if(!is.null(AIC.list)){
    bestmodind<-which(AIC.list==min(AIC.list))
    if(length(bestmodind)>1){
      bestmodind<-sample(bestmodind,1)
    }
    bestmod<-mod.list[[bestmodind]]
    cfs<-coef(bestmod)
    expected<-nbcurve(dat$temperature,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
    rsqr<-1-sum((dat$growth.rate-expected)^2)/sum((dat$growth.rate-mean(dat$growth.rate))^2)
  }
  
  # Use the curve fit to find Topt and the estimated maximum growth rate (i.e. growth rate at Topt)
  grfunc<-function(x){
    -nbcurve(x,cfs[[1]],cfs[[2]],cfs[[3]],cfs[[4]])
  }
  optinfo<-optim(c(x=cfs[[1]]),grfunc)
  opt<-optinfo$par[[1]]
  maxgrowth<- -optinfo$value
  
  #stash results		
  rsqr.list[i]<-rsqr
  z.list[i]<-cfs[[1]]
  w.list[i]<-cfs[[2]]
  a.list[i]<-cfs[[3]]
  b.list[i]<-cfs[[4]]
  s.list[i]<-cfs[[5]]
  topt.list[i]<-opt
  maxgrowth.list[i]<-maxgrowth
  n.list[i]<-length(dat$temperature)
}


fits_2015 <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list)

#split curve.id.list into species and run columns
temp<-strsplit(fits_2015$curve.id.list, split=" ")
mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
df   <- as.data.frame(mat)
fits_2015_full   <- cbind(df, fits_2015) %>%
  rename(Species=V1, run=V2)

#for now just remove the models that hit zero, but in selection above I should be able to not allow those to be picked
fits_2015_full_no_zero <- fits_2015_full %>%
  filter(maxgrowth.list>0)

write_csv(fits_2015_full_no_zero, "data-processed/Norberg_TPC_fits_2015.csv")


fits_2015_full<-read_csv("data-processed/Norberg_TPC_fits_2015.csv")

dim(fits_2015_full)

par(mfrow=c(2,2))
for(i in 1:4){
  temppoints<-dat.full %>%
    filter(Species==unique(dat.full$Species)[i])
with(temppoints, plot(growth.rate~temperature,xlim=c(3,38), ylim=c(-2, 8), xlab='Temperature', 
     ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5))
  speciesfits<-subset(fits_2015_full, unique(fits_2015_full$Species) == i)
  tempcurves<-fits_2015_full %>%
    filter(Species==unique(dat.full$Species)[i]) 
  for(j in 1:100){
with(tempcurves, curve(nbcurve(x,z.list[j],w.list[j],a.list[j],b.list[j]),col='red', lwd=2,add=TRUE))
}
}


#build a prediction of umax for 200 temperatures between 3 and 38
newdata=seq(3, 38, length = 200)
umaxall=data.frame() #make an empty dataframe to stash results
for (i in 1:4){
  speciesmodel<-fits_2015_full %>%
    filter(Species==unique(Species)[i])
  for(j in speciesmodel$run){
    umaxpredict<-with(subset(speciesmodel, speciesmodel$run==j), nbcurve(newdata,z.list,w.list,a.list,b.list))
    dftemp<-data.frame(umaxpredict=umaxpredict, temperature=newdata, run=j, Species=speciesmodel$Species[1])
    umaxall<-rbind(umaxall, dftemp)
  }
}

head(umaxall)

write.csv(umaxall, "data-processed/Norberg_model_umaxall_predictions.csv")
#
#not sure if I need the rest of the code below
#

plot(dat.full$growth.rate[as.numeric(dat.full$curve.id) == 2]~dat.full$temperature[as.numeric(dat.full$curve.id) == 2], main=curve.id.list[2],
     xlim=c(3,38),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[2],fits$w.list[2],fits$a.list[2],fits$b.list[2]),col='red', lwd=2,add=TRUE)

plot(dat.full$growth.rate[as.numeric(dat.full$curve.id) == 3]~dat.full$temperature[as.numeric(dat.full$curve.id) == 3], main=curve.id.list[3],
     xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[3],fits$w.list[3],fits$a.list[3],fits$b.list[3]),col='red', lwd=2,add=TRUE)

plot(dat.full$growth.rate[as.numeric(dat.full$curve.id) == 4]~dat.full$temperature[as.numeric(dat.full$curve.id) == 4], main=curve.id.list[4],
     xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[4],fits$w.list[4],fits$a.list[4],fits$b.list[4]),col='red', lwd=2,add=TRUE)




