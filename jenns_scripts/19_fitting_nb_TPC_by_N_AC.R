### fit TPC to 2015 umax data 
#goal:fit Norberg model to r estimate for 1000 bootstraps 
#for each of 7 nitrate concentations for each of 7 species
#worried this will take a long time - yes it is
#therefore, fitting TPC to all the jacknkifed data - might just keep this
#or, I could use the output to set parameters for the quick fit methods of all 100 jackknifes
#quick fit = reduce the searching scope of starting conditions

#read in the data
TT_unique_boots<-read_csv("data-processed/TT_unique_boots.csv")
CH_unique_boots<-read_csv("data-processed/CH_unique_boots.csv")
CS_unique_boots<-read_csv("data-processed/CS_unique_boots.csv")
AC_unique_boots<-read_csv("data-processed/AC_unique_boots.csv")

#libraries
library(bbmle)
library(tidyverse)
library(stringr)

#plot the bootstrapped data - TPC by N
AC_unique_boots %>% 
  group_by(Temperature, N.Treatment) %>% 
  ggplot(aes(x = Temperature, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  theme_bw() + 
  stat_summary(mapping = aes(x = Temperature, y = a, group = N.Treatment),
               fun.ymin = function(z) { quantile(z,0.05) },
               fun.ymax = function(z) { quantile(z,0.95) },
               fun.y = median, pch=1, size = 0.5, geom="line")


head(AC_unique_boots)

#get data ready - just one species to start with
dat.full <- AC_unique_boots %>% 
  filter(N.Treatment!=0) %>% 
  #mutate(eachcurve=paste(N.Treatment,run)) %>%
  rename(curve.id = N.Treatment,
         temperature = Temperature,
         growth.rate = a)

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
  avals<-seq(-0.2,1.2,0.1)
  bvals<-seq(-0.2,0.3,0.1)
  
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


fits_each_N_AC <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list)

#add species and run columns - do this carefully and check
#df   <- dat.full[1:700,3:4]
#fits_each_N_CS_full   <- cbind(df, fits_each_N_CS)
write_csv(fits_each_N_AC, "data-processed/Norberg_fits_by_N_AC.csv")
head(dat.full)


#for now just remove the models that run along zero, but in selection above I should be able to not allow those to be picked
fits_each_AC_full_no_zero <- fits_each_N_AC %>%
  filter(maxgrowth.list>0.01)

par(mfrow=c(3,3))
for(i in c(1,3,5,7,2,4,6)){
  temppoints<-dat.full %>%
    filter(curve.id==unique(dat.full$curve.id)[i])
  with(temppoints, plot(growth.rate~temperature,xlim=c(3,38), ylim=c(-2, 8), xlab='Temperature', 
                        ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5))
  tempcurves<-fits_each_AC_full_no_zero %>%
    filter(curve.id.list==unique(dat.full$curve.id)[i]) 
  for(j in 1:100){
    with(tempcurves, curve(nbcurve(x,z.list[j],w.list[j],a.list[j],b.list[j]),col='red', lwd=2,add=TRUE))
  }
}

#build a prediction of umax for 200 temperatures between 3 and 38
newdata=seq(3, 38, length = 200)
r_pred_all=data.frame() #make an empty dataframe to stash results
for (i in c(1,3,5,7,2,4,6)){
  nitratemodel<-fits_each_AC_full_no_zero %>%
    filter(curve.id.list==unique(curve.id.list)[i])
  r_predict<-with(nitratemodel, nbcurve(newdata,z.list,w.list,a.list,b.list))
  dftemp<-data.frame(r_predict=r_predict, temperature=newdata, N.Treatment=nitratemodel$curve.id.list[1])
  r_pred_all<-rbind(r_pred_all, dftemp)
}


head(r_pred_all)

write.csv(r_pred_all, "data-processed/Norberg_fits_by_N_AC_predictions.csv")
