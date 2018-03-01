library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

TT <- read_csv("data/monod_data/TT_15_08_24.csv")


#take raw data, add actual N concentrations
TTN<- TT %>% 	
  filter(N.Treatment!=5 | Temperature!=28 | MINUTE!=22) %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
  mutate(N.Treatment = str_replace(N.Treatment, "1", "11")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "2", "22")) %>% 
  mutate(N.Treatment = str_replace(N.Treatment, "3", "33")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "4", "55")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "5", "110")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "6", "220")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "7", "330")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "8", "440")) %>%  
  mutate(N.Treatment = str_replace(N.Treatment, "1105", "55")) %>% 
  mutate(N.Treatment = as.numeric(N.Treatment))

#filter out data where K seems to have been reached beyond the variability expected by noise (quite arbitrary)
TTfiltered<- TT %>% 
	mutate(day = Hours.since.Innoc/24) %>% 
  filter(N.Treatment!=5 | Temperature!=28 | MINUTE!=22) %>% 
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
	filter( Temperature == 13 & N.Treatment == 8 |
	          Temperature == 13 & N.Treatment == 7 |
	          Temperature == 13 & N.Treatment == 6 | 
	          Temperature == 13 & N.Treatment == 5 |
	          Temperature == 13 & N.Treatment == 4 |
	          Temperature == 13 & N.Treatment == 3 |
	          Temperature == 13 & N.Treatment == 2 |
	          Temperature == 13 & N.Treatment == 1 & day<6 |
	          Temperature == 13 & N.Treatment == 0 & day<6 |
	        Temperature == 16 & N.Treatment == 8 |
				 	Temperature == 16 & N.Treatment == 7 |
				 	Temperature == 16 & N.Treatment == 6 | 
				 	Temperature == 16 & N.Treatment == 5  |
				 	Temperature == 16 & N.Treatment == 4 |
				 	Temperature == 16 & N.Treatment == 3 |
				 	Temperature == 16 & N.Treatment == 2 |
				 	Temperature == 16 & N.Treatment == 1 & day<6 |
				 	Temperature == 16 & N.Treatment == 0 & day<6|
						Temperature == 19 & N.Treatment == 8  |
						Temperature == 19 & N.Treatment == 7  |
						Temperature == 19 & N.Treatment == 6 | 
						Temperature == 19 & N.Treatment == 5 |
						Temperature == 19 & N.Treatment == 4 |
						Temperature == 19 & N.Treatment == 3 |
						Temperature == 19 & N.Treatment == 2 |
						Temperature == 19 & N.Treatment == 1 |
						Temperature == 19 & N.Treatment == 0 |
					Temperature == 22 & N.Treatment == 8 | 
					Temperature == 22 & N.Treatment == 7 |
					Temperature == 22 & N.Treatment == 6 | 
					Temperature == 22 & N.Treatment == 5 |
					Temperature == 22 & N.Treatment == 4 |
					Temperature == 22 & N.Treatment == 3 & day<4.5|
					Temperature == 22 & N.Treatment == 2 |
					Temperature == 22 & N.Treatment == 1 |
					Temperature == 22 & N.Treatment == 0 | 		
						  Temperature == 25 & N.Treatment == 8 |
						  Temperature == 25 & N.Treatment == 7 |
						  Temperature == 25 & N.Treatment == 6 | 
						  Temperature == 25 & N.Treatment == 5 |
						  Temperature == 25 & N.Treatment == 4 & day<3.5 |
						  Temperature == 25 & N.Treatment == 3 & day<3.5 |
						  Temperature == 25 & N.Treatment == 2 |
						  Temperature == 25 & N.Treatment == 1 & day<3.5|
						  Temperature == 25 & N.Treatment == 0 & day<3.5| 		
								Temperature == 28 & N.Treatment == 8 |
								Temperature == 28 & N.Treatment == 7 |
								Temperature == 28 & N.Treatment == 6 | 
								Temperature == 28 & N.Treatment == 5 |
								Temperature == 28 & N.Treatment == 4 & day<3 |
								Temperature == 28 & N.Treatment == 3 & day<3 |
								Temperature == 28 & N.Treatment == 2 |
								Temperature == 28 & N.Treatment == 1 & day<4|
								Temperature == 28 & N.Treatment == 0 & day<3.8) 


#take filtered data, add actual N concentrations, plot monod curves
TTfilteredN<- TTfiltered %>% 	
	mutate(N.Treatment = str_replace(N.Treatment, "1", "11")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "2", "22")) %>% 
	mutate(N.Treatment = str_replace(N.Treatment, "3", "33")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "4", "55")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "5", "110")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "6", "220")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "7", "330")) %>%
	mutate(N.Treatment = str_replace(N.Treatment, "8", "440")) %>%  
	mutate(N.Treatment = str_replace(N.Treatment, "1105", "55")) %>% 
	mutate(N.Treatment = as.numeric(N.Treatment))

write_csv(TTfilteredN,"data-processed/TTfilteredN.csv")

#fit a linear model
#fit just lm
TTfilteredN %>%
  ggplot(aes(y=log.Particles.per.ml, x=day))  +
  facet_grid(N.Treatment~Temperature) +
  stat_smooth(method=lm) +
  geom_point(data=TTN, color="red", alpha=0.5) + geom_point()

lmfits<-TTfilteredN %>%
  group_by(N.Treatment, Temperature) %>%
  do(tidy(lm(log.Particles.per.ml~day, data=.)))


#fit monod curve
TT_ks_umax<-lmfits %>%
  filter(term=="day")  %>%
  #filter(N.Treatment!=0) %>%
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(tidy(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
              control = nls.control(maxiter=500, minFactor=1/204800000))))


TT_ks_umax_pred<-lmfits %>%
  filter(term=="day")  %>%
  #filter(N.Treatment!=0) %>%
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(augment(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
              control = nls.control(maxiter=500, minFactor=1/204800000))))

#plot monod curve
lmfits %>%
  filter(term=="day")  %>%
  ggplot(aes(x=N.Treatment, y=estimate, color=as.factor(Temperature))) + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) + theme_bw() +
  geom_line(data=TT_ks_umax_pred, aes(x=N.Treatment, y=.fitted)) +
  facet_wrap(~Temperature)


TT_ks_umax %>%
  filter(term=="umax") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

TT_ks_umax %>%
  filter(term=="ks") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 4) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

lmfits %>%
  filter(term=="day")  %>%
  ggplot(aes(x=Temperature, y=estimate, color=as.factor(N.Treatment))) + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) + theme_bw() +
  facet_wrap(~N.Treatment)

#pick any 2 points and resample
TTfilteredN <- TTfilteredN %>% 	
  arrange(N.Treatment, Temperature, day)

lines<-data.frame()
for(i in 1:100){
  for (j in 1: length(unique(TTfilteredN$N.Treatment))){
    for (k in 1:length(unique(TTfilteredN$Temperature))){
      onefit<-TTfilteredN %>%
        filter(N.Treatment==unique(TTfilteredN$N.Treatment)[j] & 
                 Temperature==unique(TTfilteredN$Temperature)[k]) %>%
        sample_n(., 2, replace=FALSE) %>%
        do(tidy(lm(log.Particles.per.ml~day, data=.)))
      temp<-data.frame(slope=onefit[2,2],
                       intercept=onefit[1,2],
                       N.Treatment=unique(TTfilteredN$N.Treatment)[j],
                       Temperature=unique(TTfilteredN$Temperature)[k])
      lines<-rbind(lines, temp)
    }
  }
}

grouped_lines<-lines %>%
  group_by(N.Treatment, Temperature) %>%
  summarize(median_slope=median(slope), median_intercept=median(intercept)) 

TTfilteredN %>%
  ggplot(aes(y=log.Particles.per.ml, x=day)) + geom_point() +
  geom_abline(data=grouped_lines, aes(slope=median_slope, intercept=median_intercept)) +
  facet_grid(N.Treatment~Temperature) +
  stat_smooth(method=lm)


#try fitting exponential curve to each leg of data collection
#subset to only data > 3000
#fit exponential to every point above 3000, take a mean and sd of that estimate.
#this is a pseudoreplicate, but allows some estimation of sampling error.

init<-TTfilteredN %>% 	
  filter(DAY=="0") 

TTfilteredN <- TTfilteredN %>% 	
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
  arrange(N.Treatment, Temperature, day)  %>%
  group_by(N.Treatment, Temperature) %>% 
  mutate(start=head(Particles.per.ml, 1), log.start=head(log.Particles.per.ml, 1))

logmaxcells<-log(2000)
TTfittedlineperobs<-TTfilteredN %>%
  filter(log.Particles.per.ml>logmaxcells) %>% 	
  group_by(N.Treatment, Temperature, day) %>% 
  do(tidy(lm(I(log.Particles.per.ml-log.start) ~ 0 + day, data= .)))

TTfittedlineinit<-TTfilteredN %>%
  filter(log.Particles.per.ml>logmaxcells) %>% 	
  select(log.start)

TTfittedlineperobs_w_init<-cbind(TTfittedlineperobs, TTfittedlineinit)
  
TTfilteredN %>%
  ggplot(aes(x=day, y=log.Particles.per.ml)) + geom_point() + 
  facet_wrap(N.Treatment~Temperature) + 
  geom_abline(data=TTfittedlineperobs_w_init, aes(slope=estimate, intercept=log.start))

write_csv(TTfittedlineperobs_w_init, "data-processed/TT_fitted_r_per_obs.csv")


#
#
#
# working in non-logged data, fitting exponential curve

TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 400 * exp(a*day), #N0e^rt
							data= .,  start=list(a=0.7),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ungroup() %>% 
  filter(term=="a") %>%
	ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
	geom_line() + theme_bw() + facet_wrap( ~ Temperature)
ggsave("figures/TT_monod_curves.png")

TT_r <- TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(day),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) 

write_csv(TT_r, "data-processed/fitted_r_TT_from_2015.csv")


TTfilteredN %>% 
	group_by(Temperature, N.Treatment) %>% 
	do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(day),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	ggplot(aes(x = Temperature, y = estimate, color = factor(N.Treatment))) + geom_point(size = 2) +
	geom_line() + theme_bw()
ggsave("figures/TT_TPC_by_nitrate_curves.png", length=5, width=5)

#set up the bootstrapping
for (j in c(13, 16, 19, 22, 25, 28)){
  for (k in c(0, 11, 22, 33, 55, 110, 220, 440)){
    test<-subset(TTfilteredN, TTfilteredN$Temperature==j & TTfilteredN$N.Treatment==k)
    boot_r<-1:100
    for(i in 1:100){
      mod<-sample_n(test, length(test$day)-1) %>%
        do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(day),
                    data= .,  start=list(a=0.01),
                    control = nls.control(maxiter=100, minFactor=1/204800000))))
      boot_r[i]<-mod$estimate
      name<-data.frame(a=boot_r, Temperature=test$Temperature[1], N.Treatment=test$N.Treatment[1], run=1:100)
      assign(paste("TTboot_r_unique", test$Temperature[1], test$N.Treatment[1], sep ="_"), name)
    }
  }
}

# bind all of the objects that start with boot_r_unique_!
TT_unique_boots<-do.call("rbind", mget(ls(pattern="TTboot_r_unique"))) %>%
  select(a, Temperature, N.Treatment, run) %>%
  mutate(Species="TT")

write_csv(as.data.frame(TT_unique_boots), "data-processed/TT_unique_boots.csv")


#plot the bootstrapped data - monod
TT_unique_boots %>% 
  group_by(Temperature, N.Treatment) %>% 
  ggplot(aes(x = N.Treatment, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature) + 
  stat_summary(mapping = aes(x = N.Treatment, y = a),
    fun.ymin = function(z) { quantile(z,0.05) },
    fun.ymax = function(z) { quantile(z,0.95) },
    fun.y = median, pch=1, size = 0.5, colour="black")
ggsave("figures/TT_monod_bootstrap.png")

#plot the bootstrapped data - TPC by N
TT_unique_boots %>% 
  group_by(Temperature, N.Treatment) %>% 
  ggplot(aes(x = Temperature, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  theme_bw() + 
  stat_summary(mapping = aes(x = Temperature, y = a, group = N.Treatment),
               fun.ymin = function(z) { quantile(z,0.05) },
               fun.ymax = function(z) { quantile(z,0.95) },
               fun.y = median, pch=1, size = 0.5, geom="line")

#fit monod curve to each bootstrap
TT_ks_umax_boot<-TT_unique_boots %>%
  filter(N.Treatment!=0) %>%
  group_by(Temperature, run) %>% 
  mutate(r_estimate=a) %>% 
  do(tidy(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0, d=0),
              control = nls.control(maxiter=500, minFactor=1/204800000))))

write_csv(TT_ks_umax_boot, "data-processed/TT_ks_umax_boot.csv")


#plot monod curves with fitted line
TT_ks_umax_fitted<-TT_unique_boots %>%
  filter(N.Treatment!=0) %>%
  group_by(Temperature, run) %>% 
  mutate(r_estimate=a) %>% 
  do(augment(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
                 data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0, d=0),
                 control = nls.control(maxiter=100, minFactor=1/204800000))))

TT_unique_boots %>% 
  group_by(Temperature) %>% 
  ggplot(aes(x = N.Treatment, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  theme_bw() + facet_wrap( ~ Temperature) + 
  geom_line(data=TT_ks_umax_fitted, aes(x=N.Treatment, y=.fitted, color=as.factor(run))) 
ggsave("figures/TT_monod_fitted_bootstrapped.png")

#Joey - can you help make this nice?

#get range of ks and umax for each temp
TT_summ_ks_umax<-TT_ks_umax_boot %>%
  #filter(term=="ks") %>% 
  group_by(Temperature, term) %>% 
  summarize(., mn=median(estimate), uci=quantile(estimate, 0.95), lci=quantile(estimate, 0.05))

TT_summ_ks_umax %>%
  ggplot(aes(x = Temperature, y = mn)) + geom_point(size = 2) +
  facet_wrap( ~ term, scales = "free") +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.2)
ggsave("figures/TT_ks_umax_boot.png")


#previous coding working in linear data

#plot monod curves with error
linear_r<-TTfilteredN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+100)) %>% 
  filter(day>0.5) %>% 
  group_by(N.Treatment, Temperature) %>% 
  do(tidy(lm(Particles.per.ml ~ day, data=.))) 
linear_r %>% 
  filter(term=="day") %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  theme_bw() + facet_wrap( ~ Temperature) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

write_csv(linear_r, "data-processed/TT_r_lm.csv")
ggsave("figures/TT_monod_curves.png")

#plot linear fit over raw data
TTlinear_r_aug<- TTfilteredN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  #filter(day>1.5) %>% 
  group_by(N.Treatment, Temperature) %>% 
  do(augment(lm(Particles.per.ml ~ day, data=.))) 

TTN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  #filter(N.Treatment==55) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + 
  geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + 
  geom_line(data=subset(TTlinear_r_aug, TTlinear_r_aug$N.Treatment==55), aes(x=day, y=.fitted, colour=factor(N.Treatment)))
#geom_line(data=linear_r_aug, aes(x=day, y=.fitted, colour=factor(N.Treatment)))

TTN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  filter(Temperature==28) %>% 
  ggplot(data = ., aes(x = DAY, y = Particles.per.ml, color = factor(N.Treatment))) + 
  geom_point() +
  facet_wrap( ~ N.Treatment) + geom_line() + 
  geom_line(data=subset(TTlinear_r_aug, TTlinear_r_aug$Temperature==28), aes(x=day, y=.fitted, colour=factor(N.Treatment)))
#geom_line(data=linear_r_aug, aes(x=day, y=.fitted, colour=factor(N.Treatment)))

ggsave("figures/TT_linear_R_fits.png")


#fitting monod curve -------
TT_r_lm<-read_csv("data-processed/TT_r_lm.csv") #read in the data

TT_ks_umax<-TT_r_lm %>%
  filter(N.Treatment!=0) %>%
  filter(term=="day") %>% 
  #filter(Temperature==22) %>% 
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(tidy(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0, d=0),
              control = nls.control(maxiter=500, minFactor=1/204800000))))
write_csv(TT_ks_umax, "data-processed/TT_ks_umax.csv")


TT_ks_umax_fitted<-TT_r_lm %>%
  filter(N.Treatment!=0) %>%
  filter(term=="day") %>% 
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(augment(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
                 data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0, d=0),
                 control = nls.control(maxiter=100, minFactor=1/204800000))))

#plot monod curves with fitted line
linear_r %>% 
  filter(term=="day") %>% 
  filter(N.Treatment!=0) %>%
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  theme_bw() + facet_wrap( ~ Temperature) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) + 
  geom_line(data=TT_ks_umax_fitted, aes(x=N.Treatment, y=.fitted), color=1) 
ggsave("figures/TT_monod_fitted.png")

#stat_function(fun = function(x) test$umax[3]*(x / (test$ks[3]+ x)))
#joey please help me draw line from a dataframe that parse out over plots by temperature

TT_ks_umax %>%
  filter(term=="umax") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

TT_ks_umax %>%
  filter(term=="ks") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 4) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)


# plot R-star ------------
#R* = mKs / umax - m
#
#set m at 0.1
TT_umax<- TT_ks_umax %>%
  filter(term=="umax") 
TT_ks<- TT_ks_umax %>%
  filter(term=="ks") 
m<-rep(0.1, 6)
TT_umax$rstar<-m * TT_ks$estimate / (TT_umax$estimate - m)
TT_umax %>%
  ggplot(aes(x = Temperature, y = rstar)) + geom_point(size = 4) 

#next redo this by pulling 1000 iterations from umax and k estimates with norm dist with se
#write a function to draw random numbers with a boundary of 0:
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

TT_rstar_norm_summ<-data.frame(median=1:6, lci=1:6, uci=1:6, Temperature=TT_ks$Temperature)
for(i in 1:length(TT_umax$estimate)){
  rstar_norm<-(rtnorm(n=100, mean=TT_ks$estimate[i], sd=TT_ks$std.error[i], a=0, b=Inf)*m)/
    (rnorm(100, mean=TT_umax$estimate[i], sd=TT_umax$std.error[i])-m)
  TT_rstar_norm_summ[i,c(1:3)]<-quantile(rstar_norm, c(0.5, 0.05, 0.95))
}
write_csv(TT_rstar_norm_summ, "data-processed/TT_rstar.csv")
