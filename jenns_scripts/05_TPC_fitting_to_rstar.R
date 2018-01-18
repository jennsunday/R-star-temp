### fit TPC to R-star growth rate data


library(bbmle)
library(tidyverse)

## load data

cell_results<-read_csv("data-processed/logistic_growth_fits_r-star.csv")

dat.full <- cell_results %>% 
  mutate(corrected_r=r+0.1) %>% 
  filter(species !=3 | temp != 3) %>% 
  filter(temp != 38) %>% 
	rename(curve.id = species,
				 growth.rate = corrected_r, 
				 temperature = temp)


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
	
	## This loop fits the model using a range of different starting guesses. We choose the best one using AIC. This helps find good solutions even if there are
	# convergence problems.
	# Starting estimates for parameters 'a' and 'b' use a plausible range but with broadly spaced estimates to speed up fitting. 
	avals<-seq(-0.2,1.2,0.1)		
	bvals<-seq(-0.2,0.3,0.05)
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
	
	# If the quick fit yielded poor results (low R^2), try a more thorough search through parameter space
	if(rsqr<0.95){
		avals<-seq(-0.2,1.2,0.02)
		bvals<-seq(-0.2,0.3,0.02)
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

fits <-data.frame(curve.id.list, topt.list,maxgrowth.list,z.list,w.list,a.list,b.list,rsqr.list,s.list,n.list)
write_csv(fits, "data-processed/TPC_fits.csv")


### plot the curves
fits <- read_csv("data-processed/TPC_fits.csv")

par(mfrow=c(2,2))
plot(dat.full$growth.rate[dat.full$curve.id == 1]~dat.full$temperature[dat.full$curve.id == 1],ylim=c(-0.5,2), main=curve.id.list[1],
		 xlim=c(1,34),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[1],fits$w.list[1],fits$a.list[1],fits$b.list[1]),col='red', lwd=2,add=TRUE)

plot(dat.full$growth.rate[dat.full$curve.id == 2]~dat.full$temperature[dat.full$curve.id == 2],ylim=c(-0.5,3),main=curve.id.list[2],
		 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[2],fits$w.list[2],fits$a.list[2],fits$b.list[2]),col='red', lwd=2,add=TRUE)

plot(dat.full$growth.rate[dat.full$curve.id == 3]~dat.full$temperature[dat.full$curve.id == 3],ylim=c(-1.5,2), main=curve.id.list[3],
		 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[4],fits$w.list[4],fits$a.list[4],fits$b.list[4]),col='red', lwd=2,add=TRUE)


plot(dat.full$growth.rate[dat.full$curve.id == 4]~dat.full$temperature[dat.full$curve.id == 4],ylim=c(-0.5,2), main=curve.id.list[4],
		 xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[3],fits$w.list[3],fits$a.list[3],fits$b.list[3]),col='red', lwd=2,add=TRUE)

par(mfrow=c(1,1))
par(mar=c(4,5,0.1, 0.1))
plot(dat.full$growth.rate[dat.full$curve.id == 1]~dat.full$temperature[dat.full$curve.id == 1],ylim=c(-1.5,3), col=1, las=1, 
     xlim=c(3,38),xlab='Temperature',ylab='Intrinsic growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[1],fits$w.list[1],fits$a.list[1],fits$b.list[1]),col=1, lwd=2,add=TRUE)

points(dat.full$growth.rate[dat.full$curve.id == 2]~dat.full$temperature[dat.full$curve.id == 2],ylim=c(-0.5,3),
     xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature', col=2, ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[2],fits$w.list[2],fits$a.list[2],fits$b.list[2]), col=2, lwd=2,add=TRUE)

points(dat.full$growth.rate[dat.full$curve.id == 3]~dat.full$temperature[dat.full$curve.id == 3],ylim=c(-1.5,2), 
     xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature', col=3, ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[4],fits$w.list[4],fits$a.list[4],fits$b.list[4]), col=3, lwd=2,add=TRUE)
###hmmm... why does the fitted curve from species 4 fit the data from species 3 and vice-versa?
##Joey any thoughts?
points(dat.full$growth.rate[dat.full$curve.id == 4]~dat.full$temperature[dat.full$curve.id == 4],ylim=c(-0.5,2),
     xlim=c(min(dat$temperature)-2,max(dat$temperature)+2),xlab='Temperature',col=4, ylab='Specific growth rate (per day)', cex.lab=1.5,cex.axis=1.5)
curve(nbcurve(x,fits$z.list[3],fits$w.list[3],fits$a.list[3],fits$b.list[3]),col=4, lwd=2,add=TRUE)
abline(0,0)


