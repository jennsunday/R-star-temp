#libraries
library(bbmle)
library(tidyverse)
library(stringr)


r_pred_all_CS<-cbind(read_csv("data-processed/Norberg_fits_by_N_CS_predictions.csv"), Species="CS")
r_pred_all_AC<-cbind(read_csv("data-processed/Norberg_fits_by_N_AC_predictions.csv"), Species="AC")
r_pred_all_TT<-cbind(read_csv("data-processed/Norberg_fits_by_N_TT_predictions.csv"), Species="TT")
r_pred_all_CH<-cbind(read_csv("data-processed/Norberg_fits_by_N_CH_predictions.csv"), Species="CH")
r_pred_all<-rbind(r_pred_all_CS, r_pred_all_AC, r_pred_all_TT, r_pred_all_CH)

curvefits_CS<-cbind(read_csv("data-processed/Norberg_fits_by_N_CS.csv"), Species="CS")
curvefits_AC<-cbind(read_csv("data-processed/Norberg_fits_by_N_AC.csv"), Species="AC")
curvefits_TT<-cbind(read_csv("data-processed/Norberg_fits_by_N_TT.csv"), Species="TT")
curvefits_CH<-cbind(read_csv("data-processed/Norberg_fits_by_N_CH.csv"), Species="CH")
TPC_fits_all<-rbind(curvefits_CS, curvefits_AC, curvefits_TT, curvefits_CH) %>%
  rename(N.Treatment=curve.id.list)

head(TPC_fits_all)
head(r_pred_all)

r_pred_all %>%
  ggplot(aes(y=r_predict, x=temperature, color=as.factor(N.Treatment))) + 
  coord_cartesian(ylim = c(0, 8)) + geom_line(size=1) + facet_grid(~Species) +
  theme_bw() + 
  geom_point(data=TPC_fits_all, aes(x=topt.list, y=maxgrowth.list, fill=as.factor(N.Treatment)), 
             colour="black", pch=21)
  ggsave("figures/TPC_by_N.pdf", width=7, height=2)
  #geom_point(data=TPC_fits_all, aes(x=topt.list, y=maxgrowth.list)) +

TPC_fits_all %>%
  ggplot(aes(y=topt.list, x=N.Treatment)) + facet_grid(~Species) + geom_point()

r_pred_all_CS %>%
  ggplot(aes(y=r_predict, x=temperature, color=as.factor(N.Treatment))) + 
  coord_cartesian(ylim = c(0, 8)) + geom_line(size=1) + 
  theme_bw()


head(r_pred_all_AC)
  group_by(N.Treatment, temperature) %>%
  ggplot(aes(y=r_predict, x=temperature, color=as.factor(N.Treatment))) + 
  coord_cartesian(ylim = c(0, 8)) + geom_line(size=1) + 
  theme_bw()


r_pred_all_TT %>%
  ggplot(aes(y=r_predict, x=temperature, color=as.factor(N.Treatment))) + 
  coord_cartesian(ylim = c(0, 8)) + geom_line(size=1) + 
  theme_bw()

r_pred_all_CH %>%
  ggplot(aes(y=r_predict, x=temperature, color=as.factor(N.Treatment))) + 
  coord_cartesian(ylim = c(0, 8)) + geom_line(size=1) + 
  theme_bw()


####old code for run of TT when I iterated each jackknife