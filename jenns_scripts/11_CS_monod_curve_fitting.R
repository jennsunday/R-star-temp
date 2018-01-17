library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)
library(dplyr)

CS <- read_csv("data/monod_data/CS_15_08_04.csv")


#take raw data, add actual N concentrations
CSN<- CS %>% 	
  mutate(day = Hours.since.Innoc/24) %>% 
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

#plot raw data 
CSN %>% 
  #mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


#filter out data where K seems to have been reaCSed beyond the variability expected by noise (quite arbitrary)
CSfiltered<- CS %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  filter(  Temperature == 13 & N.Treatment == 8 & day<6.5|
           Temperature == 13 & N.Treatment == 7 & day<6.5|
           Temperature == 13 & N.Treatment == 6 & day<6.5 & day>2| 
           Temperature == 13 & N.Treatment == 5 & day<6.5|
           Temperature == 13 & N.Treatment == 4 & day<6.5|
           Temperature == 13 & N.Treatment == 3 & day<6.5|
           Temperature == 13 & N.Treatment == 2 & day<7 & day>2|
           Temperature == 13 & N.Treatment == 1 & day<7 & day>2|
           Temperature == 13 & N.Treatment == 0 & day<5.5 & day>0.5|
           Temperature == 16 & N.Treatment == 8 & day<6|
           Temperature == 16 & N.Treatment == 7 & day<6|
           Temperature == 16 & N.Treatment == 6 & day<6| 
           Temperature == 16 & N.Treatment == 5 & day<6|
           Temperature == 16 & N.Treatment == 4 & day<6|
           Temperature == 16 & N.Treatment == 3 & day<6|
           Temperature == 16 & N.Treatment == 2 & day<6|
           Temperature == 16 & N.Treatment == 1 & day<5.5|
           Temperature == 16 & N.Treatment == 0 & day<5.5 & day>0.5|
           Temperature == 19 & N.Treatment == 8 & day<6|
           Temperature == 19 & N.Treatment == 7 & day<6|
           Temperature == 19 & N.Treatment == 6 & day<6| 
           Temperature == 19 & N.Treatment == 5 & day<6|
           Temperature == 19 & N.Treatment == 4 & day<6|
           Temperature == 19 & N.Treatment == 3 & day<6 |
           Temperature == 19 & N.Treatment == 2 & day<6 |
           Temperature == 19 & N.Treatment == 1 & day<5 |
           Temperature == 19 & N.Treatment == 0 & day<3.5 & day>0.5|
           Temperature == 22 & N.Treatment == 8 & day<6 |
           Temperature == 22 & N.Treatment == 7 & day<6 |
           Temperature == 22 & N.Treatment == 6 & day<6 | 
           Temperature == 22 & N.Treatment == 5 & day<6 |
           Temperature == 22 & N.Treatment == 4 & day<5.5 |
           Temperature == 22 & N.Treatment == 3 & day<5.5 |
           Temperature == 22 & N.Treatment == 2 & day<5.5 |
           Temperature == 22 & N.Treatment == 1 & day<4 |
           Temperature == 22 & N.Treatment == 0 & day<3.5 & day>0.5|
           Temperature == 25 & N.Treatment == 8 & day<6 |
           Temperature == 25 & N.Treatment == 7 & day<6 |
           Temperature == 25 & N.Treatment == 6 & day<6 | 
           Temperature == 25 & N.Treatment == 5 & day<6 |
           Temperature == 25 & N.Treatment == 4 & day<5.5 |
           Temperature == 25 & N.Treatment == 3 & day<5.5 |
           Temperature == 25 & N.Treatment == 2 & day<5.5 |
           Temperature == 25 & N.Treatment == 1 & day<4.5 |
           Temperature == 25 & N.Treatment == 0 & day<2.8 & day>0.5|
           Temperature == 28 & N.Treatment == 8 & day<6 |
           Temperature == 28 & N.Treatment == 7 & day<6 |
           Temperature == 28 & N.Treatment == 6 & day<6 | 
           Temperature == 28 & N.Treatment == 5 & day<6 |
           Temperature == 28 & N.Treatment == 4 & day<6 |
           Temperature == 28 & N.Treatment == 3 & day<6 |
           Temperature == 28 & N.Treatment == 2 & day<6 |
           Temperature == 28 & N.Treatment == 1 & day<5.5 |
           Temperature == 28 & N.Treatment == 0 & day>0.5 & day<2.8  )

#take filtered data, add actual N concentrations
CSfilteredN <- CSfiltered %>% 	
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


#plot filtered data - linear model fits...
CSfilteredN %>% 
  filter(day>0.5) %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + 
  geom_smooth(method=lm, se=FALSE) 


linear_r<-CSfilteredN %>% 
    mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
    filter(day>0.5) %>% 
    group_by(N.Treatment, Temperature) %>% 
    do(tidy(lm(Particles.per.ml ~ day, data=.))) 
write_csv(linear_r, "data-processed/CS_r_lm.csv")

#plot monod curves with error
linear_r %>% 
    filter(term=="day") %>% 
    ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
      theme_bw() + facet_wrap( ~ Temperature) + 
    geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)
ggsave("figures/CS_monod_curves.png")

#plot linear fit over raw data
linear_r_aug<- CSfilteredN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  filter(day>0.5) %>% 
  group_by(N.Treatment, Temperature) %>% 
  do(augment(lm(Particles.per.ml ~ day, data=.))) 

CSN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  filter(N.Treatment==11) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + 
           geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + 
  geom_line(data=subset(linear_r_aug, linear_r_aug$N.Treatment==11), aes(x=day, y=.fitted, colour=factor(N.Treatment)))
  #geom_line(data=linear_r_aug, aes(x=day, y=.fitted, colour=factor(N.Treatment)))

ggsave("figures/CS_linear_R_fits.png")

#fitting monod curve -------
CS_r_lm<-read_csv("data-processed/CS_r_lm.csv") #read in the data

CS_ks_umax<-CS_r_lm %>%
  filter(N.Treatment!=0) %>%
  filter(term=="day") %>% 
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(tidy(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 0.001, umax = 0.01), algorithm="port", lower=list(c=0, d=0),
              control = nls.control(maxiter=100, minFactor=1/204800000))))
write_csv(CS_ks_umax, "data-processed/CS_ks_umax.csv")

CS_ks_umax_fitted<-CS_r_lm %>%
  filter(N.Treatment!=0) %>%
  filter(term=="day") %>% 
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(augment(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
                 data= .,  start=list(ks = 0.001, umax = 0.01), algorithm="port", lower=list(c=0, d=0),
                 control = nls.control(maxiter=100, minFactor=1/204800000))))

#plot monod curves with fitted line
linear_r %>% 
  filter(term=="day") %>% 
  filter(N.Treatment!=0) %>%
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  theme_bw() + facet_wrap( ~ Temperature) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) + 
  geom_line(data=CS_ks_umax_fitted, aes(x=N.Treatment, y=.fitted), color=1) 


#stat_function(fun = function(x) test$umax[3]*(x / (test$ks[3]+ x)))
#joey please help me draw line from a dataframe that parse out over plots by temperature

CS_ks_umax %>%
  filter(term=="umax") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)
  
CS_ks_umax %>%
  filter(term=="ks") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 4) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)


# plot R-star ------------
#R* = mKs / umax - m
#
#set m at 0.1
CS_umax<- CS_ks_umax %>%
  filter(term=="umax") 
CS_ks<- CS_ks_umax %>%
  filter(term=="ks") 
m<-rep(0.1, 6)
CS_umax$rstar<-m * CS_ks$estimate / (CS_umax$estimate - m)
CS_umax %>%
  ggplot(aes(x = Temperature, y = rstar)) + geom_point(size = 4) 

#next redo this by pulling 1000 iterations from umax and k estimates with norm dist with se
#write a function to draw random numbers with a boundary of 0:
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

CS_rstar_norm_summ<-data.frame(median=1:6, lci=1:6, uci=1:6, Temperature=CS_ks$Temperature)
for(i in 1:length(TT_umax$estimate)){
  rstar_norm<-(rtnorm(n=100, mean=CS_ks$estimate[i], sd=CS_ks$std.error[i], a=0, b=Inf)*m)/
    (rnorm(100, mean=CS_umax$estimate[i], sd=CS_umax$std.error[i])-m)
  CS_rstar_norm_summ[i,c(1:3)]<-quantile(rstar_norm, c(0.5, 0.05, 0.95))
}

write_csv(CS_rstar_norm_summ, "data-processed/CS_rstar.csv")

with(CS_rstar_norm_summ, plot(median~Temperature, ylim=c(0, max(uci)*2), ylab="R star, uM", las=1))
with(CS_rstar_norm_summ, segments(Temperature, lci, 
                                  Temperature, uci))

#
#
#
#
# previous coding working in non-logged data, fitting exponential curve
CSfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ungroup() %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature)
ggsave("figures/CS_monod_curves.png")
  

CSfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ungroup() %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature)
ggsave("figures/CS_monod_curves.png")

CS_r <- CSfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) 

write_csv(CS_r, "data-processed/fitted_r_CS_from_2015.csv")


CSfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ggplot(aes(x = Temperature, y = estimate, color = factor(Temperature))) + geom_point(size = 4) +
  geom_line() + theme_bw() + facet_wrap( ~ N.Treatment)
ggsave("figures/CS_TPC_by_nitrate_curves.png")


#set up the bootstrapping
for (j in c(13, 16, 19, 22, 25, 28)){
  for (k in c(0, 11, 22, 33, 55, 110, 220, 440)){
    test<-subset(CSfilteredN, CSfilteredN$Temperature==j & CSfilteredN$N.Treatment==k)
    boot_r<-1:1000
    for(i in 1:1000){
      mod<-sample_n(test, 5) %>% #take a sample of 5 of the sample days
        do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
                    data= .,  start=list(a=0.01),
                    control = nls.control(maxiter=100, minFactor=1/204800000))))
      boot_r[i]<-mod$estimate
      name<-data.frame(a=boot_r, Temperature=test$Temperature[1], N.Treatment=test$N.Treatment[1])
      assign(paste("CSboot_r_unique", test$Temperature[1], test$N.Treatment[1], sep ="_"), name)
    }
  }
}

par(mfrow=c(1,1))
plot(boot_r)
median(boot_r)
plot(quantile(boot_r, c(0.05, 0.5, 0.95)))


# bind all of the objects that start with boot_r_unique
CS_unique_boots<-do.call("rbind", mget(ls(pattern="CSboot_r_unique")))

#plot the jacknifed data - monod
CS_unique_boots %>% 
  group_by(Temperature, N.Treatment) %>% 
  ggplot(aes(x = N.Treatment, y = a, color = factor(N.Treatment))) + geom_point(size = 1) +
  theme_bw() + facet_wrap( ~ Temperature) + stat_summary(
    mapping = aes(x = N.Treatment, y = a),
    fun.ymin = function(z) { quantile(z,0.05) },
    fun.ymax = function(z) { quantile(z,0.95) },
    fun.y = median, pch=1, size = 0.5, colour="black")

ggsave("figures/CS_monod_bootstrapped.png")
