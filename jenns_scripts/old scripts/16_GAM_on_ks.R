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
  filter(term=="ks")  %>% 
  rename(ks=estimate)

# get data in order -------------------------------------------------------
head(all_ks) # here I have a ks estimate for 6 temperatures x 100 runs x 4 species
#now I want to gather data from each run and fit a model to ks, make a prediction for 200 temps between 3 and 38, save the output.

ksallpred=data.frame() #make an empty dataframe to stash results
for(i in 1:4){
Species_data <- all_ks %>% 
  filter(Species==unique(Species)[i])

for(j in 1:100){
  current_dataset<- Species_data %>%
    filter(run==j)
  mod<-gam(ks~s(Temperature, bs="cr", k=5), data=current_dataset)
  newdata=data.frame(Temperature=seq(3, 38, length = 200))
  pred<-predict(mod,newdata=newdata, se.fit=FALSE)
  prednozero<-ifelse(pred>0, pred, 0)
  with(current_dataset, points(ks~Temperature))
  lines(sort(newdata$Temperature), prednozero[order(newdata$Temperature)], col="#41ab5d")
  #fu <- pred$fit + pred$se.fit
  #fl <- pred$fit - pred$se.fit
  #polygon(c(newdata$Temperature, rev(newdata$Temperature)), c(fl, rev(fu)), border=FALSE, col="#41ab5d40")
  name<-data.frame(ks=prednozero, Temperature=newdata$Temperature, run=j, Species=current_dataset$Species[1])
  ksallpred<-rbind(ksallpred, name)
}
}

write_csv(ksallpredshort, "data-processed/GAM_ks_all.csv")

