library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

AC <- read_csv("data/monod_data/AC_15_08_17.csv")

#take raw data, add actual N concentrations
ACN<- AC %>% 	
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


#filter out data where K seems to have been reaACed beyond the variability expected by noise (quite arbitrary)
ACfiltered<- AC %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
  filter(log.Particles.per.ml>2.5) %>%
  filter(  Temperature == 13 & N.Treatment == 8 & day<8 |
           Temperature == 13 & N.Treatment == 7 & day<8|
           Temperature == 13 & N.Treatment == 6 & day<8| 
           Temperature == 13 & N.Treatment == 5 & day<8|
           Temperature == 13 & N.Treatment == 4 & day<8|
           Temperature == 13 & N.Treatment == 3 & day<8 |
           Temperature == 13 & N.Treatment == 2 & day<8|
           Temperature == 13 & N.Treatment == 1 & day<8 |
           Temperature == 13 & N.Treatment == 0 & day<8|
           Temperature == 16 & N.Treatment == 8 & day<8|
           Temperature == 16 & N.Treatment == 7 & day<8|
           Temperature == 16 & N.Treatment == 6 & day<8| 
           Temperature == 16 & N.Treatment == 5 & day<8|
           Temperature == 16 & N.Treatment == 4 & day<8|
           Temperature == 16 & N.Treatment == 3 & day<8|
           Temperature == 16 & N.Treatment == 2 & day<8|
           Temperature == 16 & N.Treatment == 1 & day<8 |
           Temperature == 16 & N.Treatment == 0 & day<8|
           Temperature == 19 & N.Treatment == 8 & day<8|
           Temperature == 19 & N.Treatment == 7 & day<8|
           Temperature == 19 & N.Treatment == 6 & day<8| 
           Temperature == 19 & N.Treatment == 5 & day<8|
           Temperature == 19 & N.Treatment == 4 & day<8 |
           Temperature == 19 & N.Treatment == 3 & day<8 |
           Temperature == 19 & N.Treatment == 2 & day<8 |
           Temperature == 19 & N.Treatment == 1 & day<8 |
           Temperature == 19 & N.Treatment == 0 & day<8|
           Temperature == 22 & N.Treatment == 8 & day<8 |
           Temperature == 22 & N.Treatment == 7 & day<8 |
           Temperature == 22 & N.Treatment == 6 & day<8 | 
           Temperature == 22 & N.Treatment == 5 & day<8 |
           Temperature == 22 & N.Treatment == 4 & day<8 |
           Temperature == 22 & N.Treatment == 3 & day<8 |
           Temperature == 22 & N.Treatment == 2 & day<8 |
           Temperature == 22 & N.Treatment == 1 & day<8 |
           Temperature == 22 & N.Treatment == 0 & day<8|
           Temperature == 25 & N.Treatment == 8 & day<8 |
           Temperature == 25 & N.Treatment == 7 & day<8 |
           Temperature == 25 & N.Treatment == 6 & day<8 | 
           Temperature == 25 & N.Treatment == 5 & day<8 |
           Temperature == 25 & N.Treatment == 4 & day<8 |
           Temperature == 25 & N.Treatment == 3 & day<8 |
           Temperature == 25 & N.Treatment == 2 & day<8 |
           Temperature == 25 & N.Treatment == 1 & day<8 |
           Temperature == 25 & N.Treatment == 0 & day<8|
           Temperature == 28 & N.Treatment == 8 & day<8 |
           Temperature == 28 & N.Treatment == 7 & day<8 |
           Temperature == 28 & N.Treatment == 6 & day<8 | 
           Temperature == 28 & N.Treatment == 5 & day<8 |
           Temperature == 28 & N.Treatment == 4 & day<8 |
           Temperature == 28 & N.Treatment == 3 & day<8 |
           Temperature == 28 & N.Treatment == 2 & day<8 |
           Temperature == 28 & N.Treatment == 1 & day<8 |
           Temperature == 28 & N.Treatment == 0 & day<8)


#take filtered data, add actual N concentrations, plot monod curves
ACfilteredN<- ACfiltered %>% 	
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

write_csv(ACfilteredN,"data-processed/ACfilteredN.csv")

#Vis linear model over raw data
ACfilteredN %>%
  ggplot(aes(y=log.Particles.per.ml, x=day))  +
  facet_grid(Temperature~N.Treatment) +
  stat_smooth(method=lm) +
  geom_point(data=ACN, color="red", alpha=0.5) + geom_point()

lmfits<-ACfilteredN %>%
  group_by(N.Treatment, Temperature) %>%
  do(tidy(lm(log.Particles.per.ml~day, data=.)))

#fit monod curve
AC_ks_umax<-lmfits %>%
  filter(term=="day")  %>%
  #filter(N.Treatment!=0) %>%
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(tidy(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
              control = nls.control(maxiter=500, minFactor=1/204800000))))


AC_ks_umax_pred<-lmfits %>%
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
  geom_line(data=AC_ks_umax_pred, aes(x=N.Treatment, y=.fitted)) +
  facet_wrap(~Temperature)


AC_ks_umax %>%
  filter(term=="umax") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

AC_ks_umax %>%
  filter(term=="ks") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 4) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)


# previous coding working in non-logged data, fitting exponential curve

ACfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(day),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ungroup() %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.1) 
ggsave("figures/AC_monod_curves.png")

AC_r <- ACfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(day),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) 

write_ACv(AC_r, "data-processed/fitted_r_AC_from_2015.ACv")


ACfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(day),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ggplot(aes(x = Temperature, y = estimate, color = factor(Temperature))) + geom_point(size = 1) +
  geom_line() + theme_bw() + facet_wrap( ~ N.Treatment) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.1) 
ggsave("figures/AC_TPC_by_nitrate_curves.png")



#set up the bootstrapping
for (j in c(13, 16, 19, 22, 25, 28)){
  for (k in c(0, 11, 22, 33, 55, 110, 220, 440)){
    test<-subset(ACfilteredN, ACfilteredN$Temperature==j & ACfilteredN$N.Treatment==k)
    boot_r<-1:100
    for(i in 1:100){
      mod<-sample_n(test, length(test$day)-1) %>% #take a sample of 1 - total # of days
        do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(day),
                    data= .,  start=list(a=0.01),
                    control = nls.control(maxiter=100, minFactor=1/204800000))))
      boot_r[i]<-mod$estimate
      name<-data.frame(a=boot_r, Temperature=test$Temperature[1], N.Treatment=test$N.Treatment[1], run=1:100)
      assign(paste("ACboot_r_unique", test$Temperature[1], test$N.Treatment[1], sep ="_"), name)
    }
  }
}

# bind all of the objects that start with boot_r_unique_!
AC_unique_boots<-do.call("rbind", mget(ls(pattern="ACboot_r_unique"))) %>%
  select(a, Temperature, N.Treatment, run) %>%
  mutate(Species="AC")

write_csv(as.data.frame(AC_unique_boots), "data-processed/AC_unique_boots.csv")

dim(AC_unique_boots)

#plot the bootstrapped data - monod
AC_unique_boots %>% 
  group_by(Temperature, N.Treatment) %>% 
  ggplot(aes(x = N.Treatment, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature) + 
  stat_summary(mapping = aes(x = N.Treatment, y = a),
               fun.ymin = function(z) { quantile(z,0.05) },
               fun.ymax = function(z) { quantile(z,0.95) },
               fun.y = median, pch=1, size = 0.5, colour="black")
ggsave("figures/AC_monod_bootstrap.png")

#fit monod curve to each bootstrap
AC_ks_umax_boot<-AC_unique_boots %>%
  filter(N.Treatment!=0) %>%
  group_by(Temperature, run) %>% 
  mutate(r_estimate=a) %>% 
  do(tidy(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0, d=0),
              control = nls.control(maxiter=500, minFactor=1/204800000))))
write_csv(AC_ks_umax_boot, "data-processed/AC_ks_umax_boot.csv")


#plot monod curves with fitted line
AC_ks_umax_fitted<-AC_unique_boots %>%
  filter(N.Treatment!=0) %>%
  group_by(Temperature, run) %>% 
  mutate(r_estimate=a) %>% 
  do(augment(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
                 data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0, d=0),
                 control = nls.control(maxiter=100, minFactor=1/204800000))))

AC_unique_boots %>% 
  group_by(Temperature) %>% 
  filter(N.Treatment!=0) %>%
  ggplot(aes(x = N.Treatment, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  theme_bw() + facet_wrap( ~ Temperature) + 
  geom_line(data=AC_ks_umax_fitted, aes(x=N.Treatment, y=.fitted, color=as.factor(run))) 
ggsave("figures/AC_monod_fitted_bootstrapped.png")
#Joey - can you help make this nice?

#get range of ks and umax for each temp
AC_summ_ks_umax<-AC_ks_umax_boot %>%
  #filter(term=="ks") %>% 
  group_by(Temperature, term) %>% 
  summarize(., mn=median(estimate), uci=quantile(estimate, 0.95), lci=quantile(estimate, 0.05))

AC_summ_ks_umax %>%
  ggplot(aes(x = Temperature, y = mn)) + geom_point(size = 2) +
  facet_wrap( ~ term, scales = "free") +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.2)
ggsave("figures/AC_ks_umax_boot.png")

#
#
#
# previous doce working in linear data...
ACfilteredN %>% 
  #filter(day>0.5) %>% 
  mutate(Particles.per.ml = log(Particles.per.ml + 1)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + 
  geom_smooth(method=lm, se=FALSE) 


linear_r<-ACfilteredN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  #filter(day>0.5) %>% 
  group_by(N.Treatment, Temperature) %>% 
  do(tidy(lm(Particles.per.ml ~ day, data=.))) 


#plot monod curves with error
linear_r %>% 
  filter(term=="day") %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  theme_bw() + facet_wrap( ~ Temperature) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

write_csv(linear_r, "data-processed/AC_r_lm.csv")
ggsave("figures/AC_monod_curves.png")

#plot linear fit over raw data
linear_r_aug<- ACfilteredN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  #filter(day>0.5) %>% 
  group_by(N.Treatment, Temperature) %>% 
  do(augment(lm(Particles.per.ml ~ day, data=.))) 

ACN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+1)) %>% 
  #filter(N.Treatment==440) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + 
  geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + 
  #geom_line(data=subset(linear_r_aug, linear_r_aug$N.Treatment==440), aes(x=day, y=.fitted, colour=factor(N.Treatment)))
  geom_line(data=linear_r_aug, aes(x=day, y=.fitted, colour=factor(N.Treatment)))

ggsave("figures/AC_linear_R_fits.png")

#fitting monod curve -------
AC_r_lm<-read_csv("data-processed/AC_r_lm.csv") #read in the data

AC_ks_umax<-AC_r_lm %>%
  filter(N.Treatment!=0) %>%
  filter(term=="day") %>% 
  group_by(Temperature) %>% 
  mutate(r_estimate=estimate) %>% 
  do(tidy(nls(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
              data= .,  start=list(ks = 0.001, umax = 0.01), algorithm="port", lower=list(c=0, d=0),
              control = nls.control(maxiter=100, minFactor=1/204800000))))
write_csv(AC_ks_umax, "data-processed/AC_ks_umax.csv")

AC_ks_umax_fitted<-AC_r_lm %>%
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
  geom_line(data=AC_ks_umax_fitted, aes(x=N.Treatment, y=.fitted), color=1) 


#stat_function(fun = function(x) test$umax[3]*(x / (test$ks[3]+ x)))
#joey please help me draw line from a dataframe that parse out over plots by temperature

AC_ks_umax %>%
  filter(term=="umax") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)

AC_ks_umax %>%
  filter(term=="ks") %>%
  ggplot(aes(x = Temperature, y = estimate)) + geom_point(size = 4) +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)


# plot R-star ------------
#R* = mKs / umax - m
#
#set m at 0.1
AC_umax<- AC_ks_umax %>%
  filter(term=="umax") 
AC_ks<- AC_ks_umax %>%
  filter(term=="ks") 
m<-rep(0.1, 6)
AC_umax$rstar<-m * AC_ks$estimate / (AC_umax$estimate - m)
AC_umax %>%
  ggplot(aes(x = Temperature, y = rstar)) + geom_point(size = 4) 

#next redo this by pulling 1000 iterations from umax and k estimates with norm dist with se
#write a function to draw random numbers with a boundary of 0:
rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}

AC_rstar_norm_summ<-data.frame(median=1:6, lci=1:6, uci=1:6, Temperature=AC_ks$Temperature)
for(i in 1:length(TT_umax$estimate)){
  rstar_norm<-(rtnorm(n=100, mean=AC_ks$estimate[i], sd=AC_ks$std.error[i], a=0, b=Inf)*m)/
    (rnorm(100, mean=AC_umax$estimate[i], sd=AC_umax$std.error[i])-m)
  AC_rstar_norm_summ[i,c(1:3)]<-quantile(rstar_norm, c(0.5, 0.05, 0.95))
}
write_csv(AC_rstar_norm_summ, "data-processed/AC_rstar.csv")

with(AC_rstar_norm_summ, plot(median~Temperature, ylim=c(0, max(uci)*1.2), ylab="R star, uM", las=1))
with(AC_rstar_norm_summ, segments(Temperature, lci, 
                                  Temperature, uci))
