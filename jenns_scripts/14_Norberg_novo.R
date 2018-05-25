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
#fitting NB to cell density
#fixed k
####################################
all_params_fit<-data.frame()
all_summary_fit<-data.frame()
all_preds_fit<-data.frame()
all_roots_ks_umax<-data.frame()

for(k in 1:4){
  data<-filter(alldata, Species==unique(alldata$Species)[k])
  fit<-nls_multstart(log.Particles.per.ml 
                        ~ N_init + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks + N.Treatment)) * day,
                        data= data,  iter = 250,
                        start_lower = c(N_init=4, a = 0.1, b=0.0001, z=5, w=5, ks = 1),
                        start_upper = c(N_init=7, a = 0.6, b=0.2, z=40, w=40, ks = 15),
                        supp_errors = 'Y',
                        convergence_count = 100,
                        na.action = na.omit,
                        lower=c(N_init=0, a = 0, b=0, z=-20, w=0, ks = 0))

  params_fit<-tidy(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_params_fit<-rbind(all_params_fit, params_fit)

  summary_fit<-glance(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_summary_fit<-rbind(all_summary_fit, summary_fit)

  preds_fit<- augment(fit)  %>% mutate(Species=unique(alldata$Species)[k])
  all_preds_fit<-rbind(all_preds_fit, preds_fit)

#get umax, ks, and roots across temp ####
#assign parameters
  a<-filter(params_fit, term=="a")$estimate
  b<-filter(params_fit, term=="b")$estimate
  z<-filter(params_fit, term=="z")$estimate
  w<-filter(params_fit, term=="w")$estimate
  ks<-filter(params_fit, term=="ks")$estimate


  #get root for every temperature
  roots_NB<-data.frame() #make an empty dataframe - temporary
  temp_NB<-data.frame() #make an empty dataframe - temporary
  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i] #for every half-degree from 1 to 50
    NB_fix<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2) * N/(ks + N) - 0.1}
    thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
    roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                           Temperature=Temperature)) #add this to the output, setting root at a large number (200) if not within interval
  }
  roots_ks_umax<-mutate(roots_NB, umax=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                ks=ks, Species=unique(alldata$Species)[k])
  all_roots_ks_umax<-rbind(all_roots_ks_umax, roots_ks_umax)
}

#View(all_roots_ks_umax)

#bootstrap this
all_params_boot<-data.frame()
roots_ks_umax_boot<-data.frame()
all_roots_ks_umax_boot<-data.frame()

for(k in 1:4){
  data<-filter(alldata, Species==unique(alldata$Species)[k])
  fit_boots <- data %>% 
    modelr::bootstrap(n = 20, id = 'boot_num') %>%
    group_by(boot_num) %>%
    mutate(fit = map(strap, ~nls_multstart(log.Particles.per.ml 
                                           ~ N_init + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks + N.Treatment)) * day,
                                           data= data.frame(.),  iter = 100,
                                           start_lower = c(N_init=4, a = 0.1, b=0.0001, z=5, w=5, ks = 1),
                                           start_upper = c(N_init=7, a = 0.6, b=0.2, z=40, w=40, ks = 15),
                                           supp_errors = 'Y',
                                           convergence_count = 100,
                                           na.action = na.omit,
                                           lower=c(N_init=0, a = 0, b=0, z=-20, w=0, ks = 0))
   ))

  # get parameters ####
  params_boot <- fit_boots %>%
   unnest(fit %>% map(tidy)) %>%
   ungroup() %>%
   mutate(Species=unique(alldata$Species)[k])

  all_params_boot<-rbind(all_params_boot, params_boot)

  #get umax, ks, and root for every bootstrap
  roots_boot<-data.frame() #make an empty dataframe
  umax_ks_boot<-data.frame() #make an empty dataframe

  for (j in unique(params_boot$boot_num)){ #for each boot
    a<-filter(params_boot, boot_num==j & term=="a")$estimate
    b<-filter(params_boot, boot_num==j & term=="b")$estimate
    z<-filter(params_boot, boot_num==j &  term=="z")$estimate
    w<-filter(params_boot, boot_num==j & term=="w")$estimate
    ks<-filter(params_boot, boot_num==j &  term=="ks")$estimate

  
    roots_NB<-data.frame() #make an empty dataframe
  
    for(i in 1:length(newdata$Temperature)){
      Temperature<-newdata$Temperature[i] #for every half-degree from 1 to 50
      NB_BR<-function(N){
         growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2) * N/(ks + N) - 0.1}
         thisroot<-try(uniroot(NB_BR, interval=c(0,200)), TRUE) #try taking the root
         roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                              Temperature=Temperature)) #add this to the output, setting root at a large number (200) if not within interval
      }
    roots_ks_umax_temp<-mutate(roots_NB, boot_num=j, 
                               umax=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                               ks=ks, Species=unique(alldata$Species)[k]) #roots, ks, and umax from one boot
    roots_ks_umax_boot<-rbind(roots_ks_umax_boot, roots_ks_umax_temp) #store all the boots together
  }
  
all_roots_ks_umax_boot<-rbind(all_roots_ks_umax_boot, roots_ks_umax_boot)
}

#View(all_roots_ks_umax_boot)

#save outputs
write_csv(all_roots_ks_umax_boot, "data-processed/all_roots_ks_umax_boot_NB_fixed.csv")
write_csv(all_roots_ks_umax, "data-processed/all_roots_ks_umax_NB_fixed.csv")
write_csv(all_summary_fit, "data-processed/all_summary_fit_NB_fixed.csv")
write_csv(all_params_fit, "data-processed/all_params_fit_NB_fixed.csv")

#read in outputs
all_roots_ks_umax_boot<-read_csv("data-processed/all_roots_ks_umax_boot_NB_fixed.csv")
all_roots_ks_umax<-read_csv("data-processed/all_roots_ks_umax_NB_fixed.csv")

#get the confidence intervals
all_boots_confint<-all_roots_ks_umax_boot %>%
  group_by(Temperature, Species) %>%
  summarise(umax_lwr_CI = quantile(umax, 0.025),
            umax_upr_CI = quantile(umax, 0.975),
            root_lwr_CI = quantile(root, 0.025),
            root_upr_CI = quantile(root, 0.975),
            ks_lwr_CI = quantile(ks, 0.025),
            ks_upr_CI = quantile(ks, 0.975)
            ) %>%
  ungroup()


write_csv(all_boots_confint, "data-processed/all_boots_confint_NB_fixed")

#plot umax
all_roots_ks_umax %>%
  ggplot() +
  coord_cartesian(ylim = c(0, 1.5), xlim = c(0, 50)) +
  theme_bw() + ylab("umax") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=umax, colour=as.factor(Species))) + 
  geom_ribbon(data=all_boots_confint, aes(x=Temperature, ymin = umax_lwr_CI, ymax = umax_upr_CI), alpha = .1)
ggsave("figures/umax_NB_fix_booted.pdf", width=8, height=2)

#plot root
all_roots_ks_umax %>%
  ggplot() + 
  coord_cartesian(ylim = c(0, 10), xlim = c(0, 50)) +
  theme_bw() + ylab("umax") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=root, colour=as.factor(Species))) +
  geom_ribbon(data=all_boots_confint, aes(x=Temperature, ymin = root_lwr_CI, ymax = root_upr_CI), alpha = .1)
ggsave("figures/root_NB_fix_booted.pdf", width=8, height=2)

all_preds_fit %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted)) + 
  geom_point() +
  facet_grid(~Species)


 ####################################
#fitting NB to cell density
#linear relationship between k and temperature
####################################
all_params_fit<-data.frame()
all_summary_fit<-data.frame()
all_preds_fit<-data.frame()
all_roots_ks_umax<-data.frame()

for(k in 1:4){
  data<-filter(alldata, Species==unique(alldata$Species)[k])
  fit<-nls_multstart(log.Particles.per.ml 
                     ~ N_init + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * 
                       (N.Treatment / (ks_a^(Temperature+ks_b) + N.Treatment)) * day,
                     data= data,  iter = 250,
                     start_lower = c(N_init=4, a = 0.1, b=0.0001, z=5, w=5, ks_a = 0, ks_b=0),
                     start_upper = c(N_init=7, a = 0.6, b=0.2, z=40, w=40, ks_a = 5, ks_b=5),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower=c(N_init=0, a = 0, b=0, z=-20, w=0, ks_a =0, ks_b=-20))
  
  params_fit<-tidy(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_params_fit<-rbind(all_params_fit, params_fit)
  
  summary_fit<-glance(fit) %>% mutate(Species=unique(alldata$Species)[k])
  all_summary_fit<-rbind(all_summary_fit, summary_fit)
  
  preds_fit<- augment(fit)  %>% mutate(Species=unique(alldata$Species)[k])
  all_preds_fit<-rbind(all_preds_fit, preds_fit)

  #get umax, ks, and roots across temp ####
  #assign parameters
  a<-filter(params_fit, term=="a")$estimate
  b<-filter(params_fit, term=="b")$estimate
  z<-filter(params_fit, term=="z")$estimate
  w<-filter(params_fit, term=="w")$estimate
  ks_a<-filter(params_fit, term=="ks_a")$estimate
  ks_b<-filter(params_fit, term=="ks_b")$estimate
  

  #get root for every temperature
  roots_NB<-data.frame() #make an empty dataframe - temporary
  for(i in newdata$Temperature){ #for temperature
    NB_fix<-function(N){
      growth_rate<-a*exp(b*i)*(1-((i-z)/(w/2))^2) * 
        N/(ks_a^(i+ks_b)  + N) - 0.1}
    thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
    roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                         Temperature=i)) #add this to the output, setting root at a large number (200) if not within interval
  }
  
  roots_ks_umax<-mutate(roots_NB, umax=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                        ks=ks_a^(newdata$Temperature+ks_b), Species=unique(alldata$Species)[k])
  all_roots_ks_umax<-rbind(all_roots_ks_umax, roots_ks_umax)
}

#View(all_roots_ks_umax)
#View(all_params_fit)

#bootstrap this
all_params_boot<-data.frame()
roots_ks_umax_boot<-data.frame()
all_roots_ks_umax_boot<-data.frame()

for(k in 1:4){
  data<-filter(alldata, Species==unique(alldata$Species)[k]) #for  one species at a time
  fit_boots <- data %>%     
    modelr::bootstrap(n = 20, id = 'boot_num') %>%           #bootstrap data n times
    group_by(boot_num) %>%
    mutate(fit = map(strap, ~nls_multstart(log.Particles.per.ml 
                                           ~N_init + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * 
                                             (N.Treatment / (ks_a^(Temperature+ks_b) + N.Treatment)) * day,
                                           data= data.frame(.),  iter = 250,
                                           start_lower = c(N_init=4, a = 0.1, b=0.0001, z=5, w=5, ks_a = 0, ks_b=-20),
                                           start_upper = c(N_init=7, a = 0.6, b=0.2, z=40, w=40, ks_a = 5, ks_b=-3),
                                           supp_errors = 'Y',
                                           convergence_count = 100,
                                           na.action = na.omit,
                                           lower=c(N_init=0, a = 0, b=0, z=-50, w=0, ks_a =0, ks_b=-50),
                                           upper=c(N_init=200, a = 200, b=200, z=200, w=200, ks_a =200, ks_b=-3))
                     ))
  
  # get parameters ####
  params_boot <- fit_boots %>%
    unnest(fit %>% map(tidy)) %>%
    ungroup() %>%
    mutate(Species=unique(alldata$Species)[k])
  
  all_params_boot<-rbind(all_params_boot, params_boot)
View(params_boot)
  #get umax, ks, and root for every bootstrap
  roots_boot<-data.frame() #make an empty dataframe
  umax_ks_boot<-data.frame() #make an empty dataframe
  
  for (j in unique(params_boot$boot_num)){ #for each boot
    a<-filter(params_boot, boot_num==j & term=="a")$estimate
    b<-filter(params_boot, boot_num==j & term=="b")$estimate
    z<-filter(params_boot, boot_num==j &  term=="z")$estimate
    w<-filter(params_boot, boot_num==j & term=="w")$estimate
    ks_a<-filter(params_boot, boot_num==j &  term=="ks_a")$estimate
    ks_b<-filter(params_boot, boot_num==j &  term=="ks_b")$estimate
    

    
    #get root for every temperature
    roots_NB<-data.frame() #make an empty dataframe - temporary
    for(i in newdata$Temperature){ #for temperatures where ks is positive
      NB_fix<-function(N){
        growth_rate<-a*exp(b*i)*(1-((i-z)/(w/2))^2) * 
          N/(ks_a^(i+ks_b)  + N) - 0.1}
      thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
      roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                           Temperature=i)) #add this to the output, setting root at a large number (200) if not within interval
    }
    


    roots_ks_umax_temp<-mutate(roots_NB, boot_num=j, 
                               umax=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                               ks=ks_a^(newdata$Temperature+ks_b), Species=unique(alldata$Species)[k]) #roots, ks, and umax from one boot
    roots_ks_umax_boot<-rbind(roots_ks_umax_boot, roots_ks_umax_temp) #store all the boots together
  }
  
  all_roots_ks_umax_boot<-rbind(all_roots_ks_umax_boot, roots_ks_umax_boot)
}

#View(all_roots_ks_umax_boot)

#save outputs
write_csv(all_roots_ks_umax_boot, "data-processed/all_roots_ks_umax_boot_NB_var.csv")
write_csv(all_roots_ks_umax, "data-processed/all_roots_ks_umax_NB_var.csv")
write_csv(all_summary_fit, "data-processed/all_summary_fit_NB_var.csv")

#read in outputs
all_roots_ks_umax_boot<-read_csv("data-processed/all_roots_ks_umax_boot_NB_var.csv")
all_roots_ks_umax<-read_csv("data-processed/all_roots_ks_umax_NB_var.csv")

#get the confidence intervals
all_boots_confint<-all_roots_ks_umax_boot %>%
  group_by(Temperature, Species) %>%
  summarise(umax_lwr_CI = quantile(umax, 0.025),
            umax_upr_CI = quantile(umax, 0.975),
            root_lwr_CI = quantile(root, 0.025),
            root_upr_CI = quantile(root, 0.975),
            ks_lwr_CI = quantile(ks, 0.025),
            ks_upr_CI = quantile(ks, 0.975)
  ) %>%
  ungroup()

write_csv(all_boots_confint, "data-processed/all_boots_confint_NB_var")


#plot umax
all_roots_ks_umax %>%
  ggplot() +
  coord_cartesian(ylim = c(0, 1.5), xlim = c(0, 50)) +
  theme_bw() + ylab("umax") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=umax, colour=as.factor(Species))) + 
  geom_ribbon(data=all_boots_confint, aes(x=Temperature, ymin = umax_lwr_CI, ymax = umax_upr_CI), alpha = .1)
ggsave("figures/umax_NB_var_booted.pdf", width=8, height=2)

#plot ks
all_roots_ks_umax %>%
  ggplot() + 
  coord_cartesian(ylim = c(0, 20), xlim = c(0, 50)) +
  theme_bw() + ylab("ks") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=ks, colour=as.factor(Species))) +
  geom_ribbon(data=all_boots_confint, aes(x=Temperature, ymin = ks_lwr_CI, ymax = ks_upr_CI), alpha = .1)
ggsave("figures/ks_NB_var_booted.pdf", width=8, height=2)

#plot root
all_roots_ks_umax %>%
  ggplot() + 
  coord_cartesian(ylim = c(0, 10), xlim = c(0, 50)) +
  theme_bw() + ylab("R-star") + facet_grid(~Species) +
  geom_line(aes(x=Temperature, y=root, colour=as.factor(Species))) +
  geom_ribbon(data=all_boots_confint, aes(x=Temperature, ymin = root_lwr_CI, ymax = root_upr_CI), alpha = .1)
ggsave("figures/root_NB_var_booted.pdf", width=8, height=2)

all_preds_fit %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted)) + 
  geom_point() +
  facet_grid(~Species)

#compare model fits
all_summary_fit_fix
all_summary_fit_var

