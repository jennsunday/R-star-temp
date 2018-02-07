#goal: fit GAM to ks data,
#to get an estiamte of ks for 200 temps between 3 and 38
#don't let ks go below zero
library(mgcv)
library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)


TT_ks_umax<-cbind(read_csv("data-processed/TT_ks_umax_boot.csv"), Species="TT")
CS_ks_umax<-cbind(read_csv("data-processed/CS_ks_umax_boot.csv"), Species="CS")
AC_ks_umax<-cbind(read_csv("data-processed/AC_ks_umax_boot.csv"), Species="AC")
CH_ks_umax<-cbind(read_csv("data-processed/CH_ks_umax_boot.csv"), Species="CH")
all_ks_umax<-rbind(TT_ks_umax, CS_ks_umax, AC_ks_umax, CH_ks_umax)
all_ks<- all_ks_umax %>%
  filter(term=="ks") 

# get data in order -------------------------------------------------------
allruns_dataset <- all_ks %>% 
  rename(ks=estimate) %>%
  filter(Species=="TT")


with(allruns_dataset, plot(ks~Temperature, type="n", xlim=c(3,38)))

for(i in 1:100){
  current_dataset<- allruns_dataset %>%
    group_by(Temperature) %>%
    sample_n(., 1)
  
  mod<-gam(ks~s(Temperature, bs="cr", k=5), data=current_dataset)
  newdata=data.frame(Temperature=seq(3, 38, length = 200))
  pred<-predict(mod,newdata=newdata, se.fit=TRUE)
  prednozero<-ifelse(pred$fit>0, pred$fit, 0)
  with(current_dataset, points(ks~Temperature))
  lines(sort(newdata$Temperature), prednozero[order(newdata$Temperature)], col="#41ab5d")
  #fu <- pred$fit + pred$se.fit
  #fl <- pred$fit - pred$se.fit
  #polygon(c(newdata$Temperature, rev(newdata$Temperature)), c(fl, rev(fu)), border=FALSE, col="#41ab5d40")
  name<-data.frame(ks=prednozero, Temperature=newdata$Temperature, run=i, Species="TT")
  assign(paste("GAM_ks", i, "TT", sep ="_"), name)
}

# bind all of the objects that start with Schoolfield_TPC
GAM_ks_all<-do.call("rbind", mget(ls(pattern="GAM_ks_")))
dim(GAM_ks_all)
head(GAM_ks_all)
write_csv(GAM_ks_all, "data-processed/GAM_ks_all.csv")

