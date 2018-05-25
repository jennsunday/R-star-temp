### fit TPC to R-star growth rate data
library(bbmle)
library(tidyverse)

## load data

cell_results<-read_csv("data-processed/logistic_growth_fits_rep_r-star.csv")%>% 
  mutate(corrected_r=r+0.1)

cell_results_less<-cell_results %>% 
  filter(species !=1 | temp != 38) %>% 
  filter(species !=2 | temp != 38) %>% 
  filter(species !=3 | temp != 38)

head(cell_results)
dat.full <- cell_results %>% 
  filter(species !=1 | temp != 38) %>% 
  filter(species !=2 | temp != 38) %>% 
  filter(species !=3 | temp != 38) %>% 
  mutate(eachcurve=paste(Species,rep)) %>%
	rename(curve.id = eachcurve,
				 growth.rate = corrected_r, 
				 temperature = temp)

#after that, rerun with all the reps combined
dat.full <- cell_results %>% 
  mutate(corrected_r=r+0.1) %>% 
  filter(species !=1 | temp != 38) %>% 
  filter(species !=2 | temp != 38) %>% 
  filter(species !=3 | temp != 38) %>% 
  mutate(eachcurve=paste(Species)) %>%
  rename(curve.id = eachcurve,
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
	avals<-seq(0.05,0.4,0.1)
	bvals<-seq(0.04,0.2,0.05)
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
	if(rsqr<0.8){
		avals<-seq(0.05,0.4,0.1)
		bvals<-seq(0.04,0.2,0.05)
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
write_csv(fits, "data-processed/TPC_fits_2017_reps_combined.csv")


### plot the curves
TPC_fits_to_2017_by_rep <- read_csv("data-processed/TPC_fits_2017.csv")
TPC_fits_to_2017 <- read_csv("data-processed/TPC_fits_2017_reps_combined.csv")

#split curve.id.list into species and run columns
#by rep
temp<-strsplit(TPC_fits_to_2017_by_rep$curve.id.list, split=" ")
mat  <- matrix(unlist(temp), ncol=2, byrow=TRUE)
df   <- as.data.frame(mat)
TPC_fits_to_2017_by_rep_full   <- cbind(df, TPC_fits_to_2017_by_rep) %>%
  rename(Species=V1, rep=V2)

#combined
temp<-strsplit(TPC_fits_to_2017$curve.id.list, split=" ")
mat  <- matrix(unlist(temp), ncol=1, byrow=TRUE)
df   <- as.data.frame(mat)
TPC_fits_to_2017_full   <- cbind(df, TPC_fits_to_2017) %>%
  rename(Species=V1)


#build a prediction of umax for 200 temperatures between 3 and 38
#by rep
newdata=seq(3, 38, length = 200)
r_pred_all=data.frame() #make an empty dataframe to stash results
for (i in 1:4){
  speciesemodel<-TPC_fits_to_2017_by_rep_full %>%
    filter(Species==unique(Species)[i])
  for(j in 1:5){
  r_predict<-with(subset(speciesemodel, speciesemodel$rep==j), nbcurve(newdata,z.list,w.list,a.list,b.list))
  dftemp<-data.frame(r_predict=r_predict, temperature=newdata, Species=speciesemodel$Species[1], rep=j)
  r_pred_all<-rbind(r_pred_all, dftemp)
}
}
write.csv(r_pred_all, "data-processed/Norberg_fits_2017_growth_by_rep.csv")

#combined
newdata=seq(3, 38, length = 200)
r_pred_all=data.frame() #make an empty dataframe to stash results
for (i in 1:4){
  speciesemodel<-TPC_fits_to_2017_full %>%
    filter(Species==unique(Species)[i])
    r_predict<-with(speciesemodel, nbcurve(newdata,z.list,w.list,a.list,b.list))
    dftemp<-data.frame(r_predict=r_predict, temperature=newdata, Species=speciesemodel$Species[1])
    r_pred_all<-rbind(r_pred_all, dftemp)
}
write.csv(r_pred_all, "data-processed/Norberg_fits_2017_growth_reps_combined.csv")


Norberg_fits_2017_by_rep<-read_csv("data-processed/Norberg_fits_2017_growth.csv")
Norberg_fits_2017_reps_combined<-read_csv("data-processed/Norberg_fits_2017_growth_reps_combined.csv")
reps_summary<-Norberg_fits_2017_by_rep %>%
  group_by(Species, temperature) %>%
  summarize(maxr=max(r_predict), minr=min(r_predict))

combined_and_min_max<-merge(Norberg_fits_2017_reps_combined, reps_summary)

combined_and_min_max %>%
  ggplot(aes(y=r_predict, x=temperature, color=as.factor(Species))) + geom_line(size=1) + 
  facet_grid(~Species) + coord_cartesian(ylim = c(0, 3)) + geom_point(data=cell_results, aes(x=temp, y=r)) +
  theme_bw()
#geom_ribbon(aes(ymax=maxr, ymin=minr, x=temperature), inherit.aes = FALSE, alpha=0.2) + 
ggsave("figures/TPC_fitted_growth_2017.pdf", width=7, height=2)

Norberg_fits_2017 %>%
  ggplot(aes(y=r_predict, x=temperature, color=Species, group=rep)) + geom_line() + 
  facet_grid(~Species) + coord_cartesian(ylim = c(0, 4)) + geom_point(data=cell_results, aes(x=temp, y=r))


cell_results_less %>%
  ggplot(aes(x=temp, y=r, color=as.factor(Species))) + geom_point() + stat_smooth(method=loess) +
  facet_grid(~Species) +  theme_bw()
ggsave("figures/loess_growth_2017.pdf", width=7, height=2)
