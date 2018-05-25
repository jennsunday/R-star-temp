#goal: get root of u(R) for every temperature in DDE + R model
#don't let ks go below zero

library(broom)
library(purrr)
library(tidyverse)
library(minpack.lm)
library(nlstools)
library(mgcv)

fit_dde<-read_csv("data-processed/fit_dde.csv")
fit_gam_pred<-read_csv("data-processed/fit_gam_pred.csv")
fit_norberg<-read_csv("data-processed/Norberg_TPC_fits_2015.csv")

newdata<-data.frame(Temperature=seq(3, 38, 0.1))

#predicting R* from IDE model with and without variation in Ks
#R-star for Species = AC
this_species<-"AC"

#parameters from IDE TPC
b1<-filter(fit_dde, Species==this_species)$estimate[1]
b2<-filter(fit_dde, Species==this_species)$estimate[2]
d0<-filter(fit_dde, Species==this_species)$estimate[3]
d1<-filter(fit_dde, Species==this_species)$estimate[4]
d2<-filter(fit_dde, Species==this_species)$estimate[5]

#parameters from Norberg TPC
z<-filter(fit_norberg, curve.id.list==this_species)$z.list
w<-filter(fit_norberg, curve.id.list==this_species)$w.list
a<-filter(fit_norberg, curve.id.list==this_species)$a.list
b<-filter(fit_norberg, curve.id.list==this_species)$b.list

#ks<-1
ks<-ifelse(filter(fit_gam_pred, Species==this_species)$.fitted>0.5, filter(fit_gam_pred, Species==this_species)$.fitted, 0.5)

#make empty data frames to stash results
roots_IDE_RB_ks_var<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks varies
roots_IDE_RB_ks_fix<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks fixed
roots_IDE_GR_ks_var<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks varies
roots_IDE_GR_ks_fix<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks fixed
roots_NB_ks_var<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks varies
roots_NB_ks_fix<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks fixed
temp_IDE_RB_ks_var<-data.frame()
temp_IDE_RB_ks_fix<-data.frame()
temp_IDE_GR_ks_var<-data.frame()
temp_IDE_GR_ks_fix<-data.frame()
temp_NB_ks_var<-data.frame()
temp_NB_ks_fix<-data.frame()


for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i]
  ks_temp<-ks[i]
  ks_fix<-mean(ks)
  
  #IDE TPC with resource limitation on birth rate
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_temp+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_var<-rbind(roots_IDE_RB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_var<-rbind(temp_IDE_RB_ks_var,Temperature)  
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_fix+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_fix<-rbind(roots_IDE_RB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_fix<-rbind(temp_IDE_RB_ks_fix,Temperature)  
  
  #IDE TPC with resource limitation on growth rate
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_temp+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_var<-rbind(roots_IDE_GR_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_var<-rbind(temp_IDE_GR_ks_var,Temperature)  
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_fix+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_fix<-rbind(roots_IDE_GR_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_fix<-rbind(temp_IDE_GR_ks_fix,Temperature)  

  #Norberg TPC with resource limitation on growth rate
  NB_R_AC_var<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_temp+N) - 0.1}
  thisroot<-try(uniroot(NB_R_AC_var, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_var<-rbind(roots_NB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_var<-rbind(temp_NB_ks_var,Temperature)  
  NB_R_AC_fix<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_fix+N) - 0.1
  }
  thisroot<-try(uniroot(NB_R_AC_fix, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_fix<-rbind(roots_NB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_fix<-rbind(temp_NB_ks_fix,Temperature)  
}

AC_roots<-data.frame(root=c(roots_IDE_RB_ks_var[,1], 
                                 roots_IDE_RB_ks_fix[,1], 
                                 roots_IDE_GR_ks_var[,1], 
                                 roots_IDE_GR_ks_fix[,1],
                                 roots_NB_ks_var[,1], 
                                 roots_NB_ks_fix[,1]),
                      temp=c(temp_IDE_RB_ks_var[,1], 
                             temp_IDE_RB_ks_fix[,1], 
                             temp_IDE_GR_ks_var[,1], 
                             temp_IDE_GR_ks_fix[,1], 
                             temp_NB_ks_var[,1], 
                             temp_NB_ks_fix[,1]),
                     model=c(rep("IDE_birth_ks_var", dim(roots_IDE_RB_ks_var)[1]),
                             rep("IDE_birth_ks_fix", dim(roots_IDE_RB_ks_fix)[1]),
                             rep("IDE_growth_ks_var", dim(roots_IDE_GR_ks_var)[1]),
                             rep("IDE_growth_ks_fix", dim(roots_IDE_GR_ks_fix)[1]),
                             rep("NB_growth_ks_var", dim(roots_NB_ks_var)[1]),
                             rep("NB_growth_ks_fix", dim(roots_NB_ks_fix)[1])),
                     Species=this_species)
  

AC_roots %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~model) +
  coord_cartesian(ylim = c(0, 100)) 
  theme_bw() + ylab("Predicted r-star, uM") 

#R-star for 
this_species<-"CH"

#parameters from IDE TPC
b1<-filter(fit_dde, Species==this_species)$estimate[1]
b2<-filter(fit_dde, Species==this_species)$estimate[2]
d0<-filter(fit_dde, Species==this_species)$estimate[3]
d1<-filter(fit_dde, Species==this_species)$estimate[4]
d2<-filter(fit_dde, Species==this_species)$estimate[5]

#parameters from Norberg TPC
z<-filter(fit_norberg, curve.id.list==this_species)$z.list
w<-filter(fit_norberg, curve.id.list==this_species)$w.list
a<-filter(fit_norberg, curve.id.list==this_species)$a.list
b<-filter(fit_norberg, curve.id.list==this_species)$b.list

#ks<-1
ks<-ifelse(filter(fit_gam_pred, Species==this_species)$.fitted>0.5, filter(fit_gam_pred, Species==this_species)$.fitted, 0.5)

#make empty data frames to stash results
roots_IDE_RB_ks_var<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks varies
roots_IDE_RB_ks_fix<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks fixed
roots_IDE_GR_ks_var<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks varies
roots_IDE_GR_ks_fix<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks fixed
roots_NB_ks_var<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks varies
roots_NB_ks_fix<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks fixed
temp_IDE_RB_ks_var<-data.frame()
temp_IDE_RB_ks_fix<-data.frame()
temp_IDE_GR_ks_var<-data.frame()
temp_IDE_GR_ks_fix<-data.frame()
temp_NB_ks_var<-data.frame()
temp_NB_ks_fix<-data.frame()


for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i]
  ks_temp<-ks[i]
  ks_fix<-mean(ks)
  
  #IDE TPC with resource limitation on birth rate
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_temp+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_var<-rbind(roots_IDE_RB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_var<-rbind(temp_IDE_RB_ks_var,Temperature)  
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_fix+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_fix<-rbind(roots_IDE_RB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_fix<-rbind(temp_IDE_RB_ks_fix,Temperature)  
  
  #IDE TPC with resource limitation on growth rate
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_temp+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_var<-rbind(roots_IDE_GR_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_var<-rbind(temp_IDE_GR_ks_var,Temperature)  
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_fix+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_fix<-rbind(roots_IDE_GR_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_fix<-rbind(temp_IDE_GR_ks_fix,Temperature)  
  
  #Norberg TPC with resource limitation on growth rate
  NB_R_AC_var<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_temp+N) - 0.1}
  thisroot<-try(uniroot(NB_R_AC_var, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_var<-rbind(roots_NB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_var<-rbind(temp_NB_ks_var,Temperature)  
  NB_R_AC_fix<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_fix+N) - 0.1
  }
  thisroot<-try(uniroot(NB_R_AC_fix, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_fix<-rbind(roots_NB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_fix<-rbind(temp_NB_ks_fix,Temperature)  
}

CH_roots<-data.frame(root=c(roots_IDE_RB_ks_var[,1], 
                            roots_IDE_RB_ks_fix[,1], 
                            roots_IDE_GR_ks_var[,1], 
                            roots_IDE_GR_ks_fix[,1],
                            roots_NB_ks_var[,1], 
                            roots_NB_ks_fix[,1]),
                     temp=c(temp_IDE_RB_ks_var[,1], 
                            temp_IDE_RB_ks_fix[,1], 
                            temp_IDE_GR_ks_var[,1], 
                            temp_IDE_GR_ks_fix[,1], 
                            temp_NB_ks_var[,1], 
                            temp_NB_ks_fix[,1]),
                     model=c(rep("IDE_birth_ks_var", dim(roots_IDE_RB_ks_var)[1]),
                             rep("IDE_birth_ks_fix", dim(roots_IDE_RB_ks_fix)[1]),
                             rep("IDE_growth_ks_var", dim(roots_IDE_GR_ks_var)[1]),
                             rep("IDE_growth_ks_fix", dim(roots_IDE_GR_ks_fix)[1]),
                             rep("NB_growth_ks_var", dim(roots_NB_ks_var)[1]),
                             rep("NB_growth_ks_fix", dim(roots_NB_ks_fix)[1])),
                     Species=this_species)


CH_roots %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~model) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() + ylab("Predicted r-star, uM") 


#Species = CS
this_species<-"CS"

#parameters from IDE TPC
b1<-filter(fit_dde, Species==this_species)$estimate[1]
b2<-filter(fit_dde, Species==this_species)$estimate[2]
d0<-filter(fit_dde, Species==this_species)$estimate[3]
d1<-filter(fit_dde, Species==this_species)$estimate[4]
d2<-filter(fit_dde, Species==this_species)$estimate[5]

#parameters from Norberg TPC
z<-filter(fit_norberg, curve.id.list==this_species)$z.list
w<-filter(fit_norberg, curve.id.list==this_species)$w.list
a<-filter(fit_norberg, curve.id.list==this_species)$a.list
b<-filter(fit_norberg, curve.id.list==this_species)$b.list

#ks<-1
ks<-ifelse(filter(fit_gam_pred, Species==this_species)$.fitted>0.5, filter(fit_gam_pred, Species==this_species)$.fitted, 0.5)

#make empty data frames to stash results
roots_IDE_RB_ks_var<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks varies
roots_IDE_RB_ks_fix<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks fixed
roots_IDE_GR_ks_var<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks varies
roots_IDE_GR_ks_fix<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks fixed
roots_NB_ks_var<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks varies
roots_NB_ks_fix<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks fixed
temp_IDE_RB_ks_var<-data.frame()
temp_IDE_RB_ks_fix<-data.frame()
temp_IDE_GR_ks_var<-data.frame()
temp_IDE_GR_ks_fix<-data.frame()
temp_NB_ks_var<-data.frame()
temp_NB_ks_fix<-data.frame()


for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i]
  ks_temp<-ks[i]
  ks_fix<-mean(ks)
  
  #IDE TPC with resource limitation on birth rate
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_temp+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_var<-rbind(roots_IDE_RB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_var<-rbind(temp_IDE_RB_ks_var,Temperature)  
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_fix+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_fix<-rbind(roots_IDE_RB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_fix<-rbind(temp_IDE_RB_ks_fix,Temperature)  
  
  #IDE TPC with resource limitation on growth rate
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_temp+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_var<-rbind(roots_IDE_GR_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_var<-rbind(temp_IDE_GR_ks_var,Temperature)  
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_fix+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_fix<-rbind(roots_IDE_GR_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_fix<-rbind(temp_IDE_GR_ks_fix,Temperature)  
  
  #Norberg TPC with resource limitation on growth rate
  NB_R_AC_var<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_temp+N) - 0.1}
  thisroot<-try(uniroot(NB_R_AC_var, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_var<-rbind(roots_NB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_var<-rbind(temp_NB_ks_var,Temperature)  
  NB_R_AC_fix<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_fix+N) - 0.1
  }
  thisroot<-try(uniroot(NB_R_AC_fix, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_fix<-rbind(roots_NB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_fix<-rbind(temp_NB_ks_fix,Temperature)  
}

CS_roots<-data.frame(root=c(roots_IDE_RB_ks_var[,1], 
                            roots_IDE_RB_ks_fix[,1], 
                            roots_IDE_GR_ks_var[,1], 
                            roots_IDE_GR_ks_fix[,1],
                            roots_NB_ks_var[,1], 
                            roots_NB_ks_fix[,1]),
                     temp=c(temp_IDE_RB_ks_var[,1], 
                            temp_IDE_RB_ks_fix[,1], 
                            temp_IDE_GR_ks_var[,1], 
                            temp_IDE_GR_ks_fix[,1], 
                            temp_NB_ks_var[,1], 
                            temp_NB_ks_fix[,1]),
                     model=c(rep("IDE_birth_ks_var", dim(roots_IDE_RB_ks_var)[1]),
                             rep("IDE_birth_ks_fix", dim(roots_IDE_RB_ks_fix)[1]),
                             rep("IDE_growth_ks_var", dim(roots_IDE_GR_ks_var)[1]),
                             rep("IDE_growth_ks_fix", dim(roots_IDE_GR_ks_fix)[1]),
                             rep("NB_growth_ks_var", dim(roots_NB_ks_var)[1]),
                             rep("NB_growth_ks_fix", dim(roots_NB_ks_fix)[1])),
                     Species=this_species)


CS_roots %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~model) +
  coord_cartesian(ylim = c(0, 100)) 
theme_bw() + ylab("Predicted r-star, uM") 

#Specoies = TT
this_species<-"TT"
#parameters from IDE TPC
b1<-filter(fit_dde, Species==this_species)$estimate[1]
b2<-filter(fit_dde, Species==this_species)$estimate[2]
d0<-filter(fit_dde, Species==this_species)$estimate[3]
d1<-filter(fit_dde, Species==this_species)$estimate[4]
d2<-filter(fit_dde, Species==this_species)$estimate[5]

#parameters from Norberg TPC
z<-filter(fit_norberg, curve.id.list==this_species)$z.list
w<-filter(fit_norberg, curve.id.list==this_species)$w.list
a<-filter(fit_norberg, curve.id.list==this_species)$a.list
b<-filter(fit_norberg, curve.id.list==this_species)$b.list

#ks<-1
ks<-ifelse(filter(fit_gam_pred, Species==this_species)$.fitted>0.5, filter(fit_gam_pred, Species==this_species)$.fitted, 0.5)

#make empty data frames to stash results
roots_IDE_RB_ks_var<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks varies
roots_IDE_RB_ks_fix<-data.frame() #store roots of IDE model with resource limitation on birth rate, ks fixed
roots_IDE_GR_ks_var<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks varies
roots_IDE_GR_ks_fix<-data.frame() #store roots of IDE model with resource limitation on growth rate, ks fixed
roots_NB_ks_var<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks varies
roots_NB_ks_fix<-data.frame() #store roots of Norberg model with resource limitation on growth rate, ks fixed
temp_IDE_RB_ks_var<-data.frame()
temp_IDE_RB_ks_fix<-data.frame()
temp_IDE_GR_ks_var<-data.frame()
temp_IDE_GR_ks_fix<-data.frame()
temp_NB_ks_var<-data.frame()
temp_NB_ks_fix<-data.frame()


for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i]
  ks_temp<-ks[i]
  ks_fix<-mean(ks)
  
  #IDE TPC with resource limitation on birth rate
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_temp+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_var<-rbind(roots_IDE_RB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_var<-rbind(temp_IDE_RB_ks_var,Temperature)  
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(ks_fix+ N))-(d0 + d1*exp(d2*Temperature)) -0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_RB_ks_fix<-rbind(roots_IDE_RB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_RB_ks_fix<-rbind(temp_IDE_RB_ks_fix,Temperature)  
  
  #IDE TPC with resource limitation on growth rate
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_temp+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_var<-rbind(roots_IDE_GR_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_var<-rbind(temp_IDE_GR_ks_var,Temperature)  
  IDE_GR<-function(N){
    growth_rate<-(b1*exp(b2*Temperature)-(d0+ d1*exp(d2*Temperature)))*(N/(ks_fix+ N)) - 0.1
  }
  thisroot<-try(uniroot(IDE_GR, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE_GR_ks_fix<-rbind(roots_IDE_GR_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE_GR_ks_fix<-rbind(temp_IDE_GR_ks_fix,Temperature)  
  
  #Norberg TPC with resource limitation on growth rate
  NB_R_AC_var<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_temp+N) - 0.1}
  thisroot<-try(uniroot(NB_R_AC_var, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_var<-rbind(roots_NB_ks_var, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_var<-rbind(temp_NB_ks_var,Temperature)  
  NB_R_AC_fix<-function(N){
    growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)*N/(ks_fix+N) - 0.1
  }
  thisroot<-try(uniroot(NB_R_AC_fix, interval=c(0,400000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_NB_ks_fix<-rbind(roots_NB_ks_fix, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_NB_ks_fix<-rbind(temp_NB_ks_fix,Temperature)  
}

TT_roots<-data.frame(root=c(roots_IDE_RB_ks_var[,1], 
                            roots_IDE_RB_ks_fix[,1], 
                            roots_IDE_GR_ks_var[,1], 
                            roots_IDE_GR_ks_fix[,1],
                            roots_NB_ks_var[,1], 
                            roots_NB_ks_fix[,1]),
                     temp=c(temp_IDE_RB_ks_var[,1], 
                            temp_IDE_RB_ks_fix[,1], 
                            temp_IDE_GR_ks_var[,1], 
                            temp_IDE_GR_ks_fix[,1], 
                            temp_NB_ks_var[,1], 
                            temp_NB_ks_fix[,1]),
                     model=c(rep("IDE_birth_ks_var", dim(roots_IDE_RB_ks_var)[1]),
                             rep("IDE_birth_ks_fix", dim(roots_IDE_RB_ks_fix)[1]),
                             rep("IDE_growth_ks_var", dim(roots_IDE_GR_ks_var)[1]),
                             rep("IDE_growth_ks_fix", dim(roots_IDE_GR_ks_fix)[1]),
                             rep("NB_growth_ks_var", dim(roots_NB_ks_var)[1]),
                             rep("NB_growth_ks_fix", dim(roots_NB_ks_fix)[1])),
                     Species=this_species)


TT_roots %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~model) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() + ylab("Predicted r-star, uM") 

allspecies_roots<-rbind(AC_roots, CH_roots, CS_roots, TT_roots)
write_csv(allspecies_roots, "data-processed/allspecies_roots.csv")

allspecies_roots %>%
  ggplot(aes(x=temp, y=root, col=model)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 100)) + 
  theme_bw() + ylab("Predicted r-star, uM")
ggsave("figures/Predicted r-star.pdf", width=7, height=2)
