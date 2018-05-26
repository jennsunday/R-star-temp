### fit TPC to 2015 umax data 
#goal:fit Norberg model to umax data for 20 bootstraps for each of 4 species


library(broom)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(tidyverse)
library(nls.multstart)

#read in data
TTfilteredN<-cbind(read_csv("data-processed/TTfilteredN.csv"), Species="TT")
CSfilteredN<-cbind(read_csv("data-processed/CSfilteredN.csv"), Species="CS")
CHfilteredN<-cbind(read_csv("data-processed/CHfilteredN.csv"), Species="CH")
ACfilteredN<-cbind(read_csv("data-processed/ACfilteredN.csv"), Species="AC")

#bring data together
alldata<-rbind(ACfilteredN, CHfilteredN, CSfilteredN, TTfilteredN)

#establish temperature frame
newdata<-data.frame(Temperature=seq(-3, 50, 0.5))

####################################
#fitting IDE to cell density
#fixed k
#ks on growth
####################################
all_params_fit<-data.frame()
all_summary_fit<-data.frame()
all_preds_fit<-data.frame()
all_roots_ks_umax<-data.frame()

for(k in 1:4){
  data<-filter(alldata, Species==unique(alldata$Species)[k])
  fit<-nls_multstart(log.Particles.per.ml 
                     ~ N_init + (((b1 * exp(b2*Temperature))
                                  -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment + zn / (ks + N.Treatment + zn)) * day,
                     data= data,  iter = 2000,
                     start_lower = c(N_init=5, b1 = 0, b2=0, d0=0, d1=0, d2=0, ks = 0, zn=0),
                     start_upper = c(N_init=7, b1 = 5, b2=2, d0=5, d1=2, d2=2, ks = 15, zn=0.5),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower=c(N_init = 0, b1 = 0, b2 = -50, d0 = 0, d1 = 0, d2 = -50, ks = 0, zn=-50),
                     control = nls.control(maxiter=1000, minFactor=1/204800000))
  
  params_fit<-tidy(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_params_fit<-rbind(all_params_fit, params_fit)

  summary_fit<-glance(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_summary_fit<-rbind(all_summary_fit, summary_fit)

  preds_fit<- augment(fit)  %>% mutate(Species=unique(alldata$Species)[k])
  all_preds_fit<-rbind(all_preds_fit, preds_fit)


#View(all_params_fit)
  #get umax, ks, and roots across temp ####
  #assign parameters
  b1<-filter(params_fit, term=="b1")$estimate
  b2<-filter(params_fit, term=="b2")$estimate
  d0<-filter(params_fit, term=="d0")$estimate
  d1<-filter(params_fit, term=="d1")$estimate
  d2<-filter(params_fit, term=="d2")$estimate
  ks<-filter(params_fit, term=="ks")$estimate
  
  #get root for every temperature
  roots_NB<-data.frame() #make an empty dataframe - temporary
  
  for(i in newdata$Temperature){#for every half-degree from 1 to 50
    NB_fix<-function(N){
      growth_rate<-((b1 * exp(b2*i)) -(d0 + d1*exp(d2*i))) * (N / (ks + N)) - 0.1}
    thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
    roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                         Temperature=i)) #add this to the output, setting root at a large number (200) if not within interval
  }
  roots_ks_umax<-mutate(roots_NB, umax=(b1 * exp(b2*newdata$Temperature))
                                         -(d0 + d1*exp(d2*newdata$Temperature)), 
                        ks=ks, Species=unique(alldata$Species)[k])
  all_roots_ks_umax<-rbind(all_roots_ks_umax, roots_ks_umax)
}

#View(all_roots_ks_umax)


#View(all_roots_ks_umax_boot)

#save outputs
write_csv(all_roots_ks_umax, "data-processed/all_roots_ks_umax_IDE_kg_fixed.csv")
write_csv(all_summary_fit, "data-processed/all_summary_fit_IDE_kg_fixed.csv")

#read in outputs
all_roots_ks_umax<-read_csv("data-processed/all_roots_ks_umax_IDE_kg_fixed.csv")


#plot umax
all_roots_ks_umax %>%
  ggplot() +
  coord_cartesian(ylim = c(0, 1.5), xlim = c(0, 50)) +
  theme_bw() + ylab("umax") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=umax, colour=as.factor(Species)))

#plot root
all_roots_ks_umax %>%
  ggplot() + 
  coord_cartesian(ylim = c(0, 10), xlim = c(0, 50)) +
  theme_bw() + ylab("umax") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=root, colour=as.factor(Species)))

all_preds_fit %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted)) + 
  geom_point() +
  facet_grid(~Species)


####################################
#fitting IDE to cell density
#ks on birth
#fixed K
####################################
all_params_fit<-data.frame()
all_summary_fit<-data.frame()
all_preds_fit<-data.frame()
all_roots_ks_umax<-data.frame()

for(k in 1:4){
  data<-filter(alldata, Species==unique(alldata$Species)[k])
  fit<-nls_multstart(log.Particles.per.ml 
                     ~ N_init + (b1 * exp(b2*Temperature) * (N.Treatment / (ks + N.Treatment))
                                  -(d0 + d1*exp(d2*Temperature))) * day,
                     data= data,  iter = 2000,
                     start_lower = c(N_init=5, b1 = 0, b2=0, ks = 0, d0=0, d1=0, d2=0),
                     start_upper = c(N_init=7, b1 = 5, b2=2, ks = 15, d0=5, d1=2, d2=2),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower = c( 0,  0, 0, 0,  0, 0, 0),
                     control = nls.control(maxiter=1000, minFactor=1/204800000))
  
  params_fit<-tidy(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_params_fit<-rbind(all_params_fit, params_fit)

  summary_fit<-glance(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_summary_fit<-rbind(all_summary_fit, summary_fit)
  
  preds_fit<- augment(fit)  %>% mutate(Species=unique(alldata$Species)[k])
  all_preds_fit<-rbind(all_preds_fit, preds_fit)
  
  #get umax, ks, and roots across temp ####
  #assign parameters
  b1<-filter(params_fit, term=="b1")$estimate
  b2<-filter(params_fit, term=="b2")$estimate
  ks<-filter(params_fit, term=="ks")$estimate
  d0<-filter(params_fit, term=="d0")$estimate
  d1<-filter(params_fit, term=="d1")$estimate
  d2<-filter(params_fit, term=="d2")$estimate

  #get root for every temperature
  roots_NB<-data.frame() #make an empty dataframe - temporary
  for(i in newdata$Temperature){ #for temperature
    NB_fix<-function(N){
      growth_rate<-((b1 * exp(b2*i))*(N / (ks + N)))-(d0 + d1 * exp(d2*i)) - 0.1}
    thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
    roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                         Temperature=i)) #add this to the output, setting root at a large number (200) if not within interval
  }
  
  roots_ks_umax<-mutate(roots_NB, umax=(b1 * exp(b2*newdata$Temperature))
                        -(d0 + d1*exp(d2*newdata$Temperature)), 
                        ks=ks, Species=unique(alldata$Species)[k])
  all_roots_ks_umax<-rbind(all_roots_ks_umax, roots_ks_umax)
}

#save outputs
write_csv(all_roots_ks_umax, "data-processed/all_roots_ks_umax_IDE_kb_fixed.csv")
write_csv(all_summary_fit, "data-processed/all_summary_fit_IDE_kb_fixed.csv")
write_csv(all_params_fit, "data-processed/all_params_fit_IDE_kb_fixed.csv")

#read in outputs
all_roots_ks_umax<-read_csv("data-processed/all_roots_ks_umax_IDE_kb_fixed.csv")

#plot umax
all_roots_ks_umax %>%
  ggplot() +
  coord_cartesian(ylim = c(0, 1.5), xlim = c(0, 50)) +
  theme_bw() + ylab("umax") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=umax, colour=as.factor(Species)))

#plot ks
all_roots_ks_umax %>%
  ggplot() + 
  coord_cartesian(ylim = c(0, 20), xlim = c(0, 50)) +
  theme_bw() + ylab("ks") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=ks, colour=as.factor(Species)))

#plot root
all_roots_ks_umax %>%
  ggplot() + 
  coord_cartesian(ylim = c(0, 10), xlim = c(0, 50)) +
  theme_bw() + ylab("R-star") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=root, colour=as.factor(Species)))

  ggsave("figures/root_NB_var_booted.pdf", width=8, height=2)

all_preds_fit %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted)) + 
  geom_point() +
  facet_grid(~Species)

#compare model fits
all_summary_fit_IDE_kb_fixed.csv<-read_csv("data-processed/all_summary_fit_IDE_kb_fixed.csv")
all_summary_fit_IDE_kg_fixed.csv<-read_csv("data-processed/all_summary_fit_IDE_kg_fixed.csv")
all_summary_fit_NB_fixed.csv<-read_csv("data-processed/all_summary_fit_NB_fixed.csv")
all_summary_fit_NB_var.csv<-read_csv("data-processed/all_summary_fit_NB_var.csv")

all_summary_fit_IDE_kb_fixed.csv
all_summary_fit_IDE_kg_fixed.csv
all_summary_fit_NB_fixed.csv
all_summary_fit_NB_var.csv
