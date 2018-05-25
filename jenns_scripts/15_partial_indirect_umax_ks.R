#Goal: use the indirect approach to predict umax ~ temp for Fig.1
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

# Model monod relationship semi-INDIRECT approach**********
#add a column of initial cell counts to each species - clunky method


N_init_fit_NB_fixed<-read_csv("data-processed/all_params_fit_NB_fixed.csv") %>%
  filter(term=="N_init")
backtogether<-data.frame()
for(i in 1:4){
subset_sp0<-N_init_fit_NB_fixed %>% 
  filter(Species==unique(alldata$Species)[i])

subset_sp<-alldata %>%
  filter(Species==unique(alldata$Species)[i]) 

subset_sp<-subset_sp%>% 
  mutate(meaninit=subset_sp0$estimate)

backtogether<-rbind(backtogether, subset_sp)
}
alldata<-backtogether

#fit model
all_params_fit<-data.frame()
all_preds_fit<-data.frame()
for(i in 1:4){
  for(j in 1:6){
  data<-alldata %>% 
        filter(Species==unique(alldata$Species)[i]) %>%
        filter(Temperature==unique(alldata$Temperature)[j])
  fit<-nls_multstart(log.Particles.per.ml 
                     ~ meaninit + (umax*(N.Treatment / (ks + N.Treatment))) * day,
                     data= data,  iter = 250,
                     start_lower = c(umax=0.1, ks = 1),
                     start_upper = c(umax=0.5, ks = 15),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower=c(umax=0, ks = 0.0001))
  
  params_fit<-tidy(fit) %>% 
    mutate(Species=unique(alldata$Species)[i]) %>% 
    mutate(Temperature=unique(alldata$Temperature)[j])
  all_params_fit<-rbind(all_params_fit, params_fit)
  
  preds_fit<- augment(fit)  %>% mutate(Species=unique(alldata$Species)[k])
  all_preds_fit<-rbind(all_preds_fit, preds_fit)
  
  }}

write_csv(all_params_fit, "data-processed/indirect_umax_ks.csv")

all_params_fit %>%
  ggplot(aes(y=estimate, x=as.numeric(Temperature))) +
    facet_grid(term~Species, scales="free_y") + geom_point() +
    geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))
         


