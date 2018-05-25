#Goal: use the indirect approach to predict growth rate for every N

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

# Model growth rate semi-INDIRECT approach**********
#set initial nitrate based on big model means
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

all_params_fit<-data.frame()
for(i in 1:4){
  for(j in 1:6){
    for(k in 1:length(unique(alldata$N.Treatment))){
    data<-alldata %>% 
      filter(Species==unique(alldata$Species)[i]) %>%
      filter(Temperature==unique(alldata$Temperature)[j])%>%
      filter(N.Treatment==unique(alldata$N.Treatment)[k])
      fit<-nls_multstart(log.Particles.per.ml 
                       ~ meaninit + r * day,
                       data= data,  iter = 250,
                       start_lower = c(r = 0.1),
                       start_upper = c(r=1.5),
                       supp_errors = 'Y',
                       convergence_count = 100,
                       na.action = na.omit,
                       lower=c(r=0))
    
    params_fit<-tidy(fit) %>% 
      mutate(Species=unique(alldata$Species)[i]) %>% 
      mutate(Temperature=unique(alldata$Temperature)[j]) %>%
      mutate(N.Treatment=unique(unique(alldata$N.Treatment)[k]))
    all_params_fit<-rbind(all_params_fit, params_fit)
  }}}

write_csv(all_params_fit, "data-processed/indirect_r.csv")


all_params_fit %>%
  ggplot(aes(y=estimate, x=as.numeric(N.Treatment), col=as.factor(Temperature))) +
  facet_grid(term~Species, scales="free_y") + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))




#get model-fitted u across temp and nitrate ####
#assign parameters

#write function
growth_function<-function (N,Temp,a,b,z,w,ks) {
  a*exp(b*Temp)*(1-((Temp-z)/(w/2))^2) * N/(ks + N)
}

predict_grid<-expand.grid(N=unique(alldata$N.Treatment), Temp=unique(alldata$Temperature), Species=unique(alldata$Species))

#add Norberg model-fitted parameters to predict_grid
NB_fit<-read_csv("data-processed/all_params_fit_NB_fixed.csv")

backtogether<-data.frame()
for(i in 1:4){
predict_grid_sp<-predict_grid %>%
  filter(Species==unique(alldata$Species)[i]) 
NB_fit_sp<-NB_fit %>%
    filter(Species==unique(alldata$Species)[i]) 
predict_grid_sp<-predict_grid_sp%>% 
    mutate(a=filter(NB_fit_sp, term=="a")$estimate,
           b=filter(NB_fit_sp, term=="b")$estimate,
           z=filter(NB_fit_sp, term=="z")$estimate,
           w=filter(NB_fit_sp, term=="w")$estimate,
           ks=filter(NB_fit_sp, term=="ks")$estimate)
  backtogether<-rbind(backtogether, predict_grid_sp)
}
predict_grid_full<-backtogether


predict_growth <- predict_grid_full %>%
  mutate(r_pred=growth_function(predict_grid_full$N,
                                predict_grid_full$Temp,
                                predict_grid_full$a,
                                predict_grid_full$b,
                                predict_grid_full$z,
                                predict_grid_full$w,
                                predict_grid_full$ks))
  

ggplot() +
  geom_line(data=predict_growth, aes(y=r_pred, x=N, colour=as.factor(Temp))) +
  facet_grid(~Species)

indirect_r_fit<-read_csv("data-processed/indirect_r.csv")
indirect_r_fit %>%
  ggplot(aes(y=estimate, x=as.numeric(N.Treatment), col=as.factor(Temperature))) +
  facet_grid(term~Species, scales="free_y") + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))+
  geom_line(data=predict_growth, aes(y=r_pred, x=N, colour=as.factor(Temp))) +
  ylab("growth rate, d-1") +
  xlab("Nitrate concentration, uM") +
  labs(col="Temperature") +
  theme_bw()
ggsave("figures/monod_fits.pdf", width=8, height=2)

indirect_r_fit %>%
  ggplot(aes(y=estimate, x=Temperature, col=as.factor(N.Treatment))) +
  facet_grid(term~Species, scales="free_y") + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))+
  geom_line(data=predict_growth, aes(y=r_pred, x=Temp, colour=as.factor(N))) +
  ylab("growth rate, d-1") +
  labs(col="Nitrate conc. uM") +
  theme_bw()
ggsave("figures/TPC_fits.pdf", width=8, height=2)

# wih IDE_kb ###########
#write function
growth_function<-function (N,Temp,b1,b2,d0,d1,d2,ks) {
  ((b1 * exp(b2*Temp))* (N / (ks + N)) -(d0 + d1*exp(d2*Temp))) 
}

predict_grid<-expand.grid(N=unique(alldata$N.Treatment), Temp=unique(alldata$Temperature), Species=unique(alldata$Species))

#add IDE_kb_fixed model-fitted parameters to predict_grid
IDE_fit<-read_csv("data-processed/all_params_fit_IDE_kb_fixed.csv")

backtogether<-data.frame()
for(i in 1:4){
  predict_grid_sp<-predict_grid %>%
    filter(Species==unique(alldata$Species)[i]) 
  IDE_fit_sp<-IDE_fit %>%
    filter(Species==unique(alldata$Species)[i]) 
  predict_grid_sp<-predict_grid_sp%>% 
    mutate(b1=filter(IDE_fit_sp, term=="b1")$estimate,
           b2=filter(IDE_fit_sp, term=="b2")$estimate,
           d0=filter(IDE_fit_sp, term=="d0")$estimate,
           d1=filter(IDE_fit_sp, term=="d1")$estimate,
           d2=filter(IDE_fit_sp, term=="d2")$estimate,
           ks=filter(IDE_fit_sp, term=="ks")$estimate)
  backtogether<-rbind(backtogether, predict_grid_sp)
}
predict_grid_full<-backtogether


predict_growth <- predict_grid_full %>%
  mutate(r_pred=growth_function(predict_grid_full$N,
                                predict_grid_full$Temp,
                                predict_grid_full$b1,
                                predict_grid_full$b2,
                                predict_grid_full$d0,
                                predict_grid_full$d1,
                                predict_grid_full$d2,
                                predict_grid_full$ks))


ggplot() +
  geom_line(data=predict_growth, aes(y=r_pred, x=N, colour=as.factor(Temp))) +
  facet_grid(~Species)


indirect_r_fit<-read_csv("data-processed/indirect_r.csv")
indirect_r_fit %>%
  ggplot(aes(y=estimate, x=as.numeric(N.Treatment), col=as.factor(Temperature))) +
  facet_grid(term~Species, scales="free_y") + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))+
  geom_line(data=predict_growth, aes(y=r_pred, x=N, colour=as.factor(Temp))) +
  theme_bw()
ggsave("figures/monod_fits.pdf", width=8, height=2)

indirect_r_fit %>%
  ggplot(aes(y=estimate, x=Temperature, col=as.factor(N.Treatment))) +
  facet_grid(term~Species, scales="free_y") + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))+
  geom_line(data=predict_growth, aes(y=r_pred, x=Temp, colour=as.factor(N))) +
  theme_bw()


