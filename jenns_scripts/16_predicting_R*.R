library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

#once I settle on a plan for curve-fitting, I need to combine both curve fits within the same 
#draw, because ks and umax are *correlated* within runs. But for now just curve fitting sorted.
#combining draws will probably reduce extremes.

#load TPC fits of umax across temperature, Gam fits of Ks across temperature, and observed Nitrate remaining in system
#TT_Schoolfield_TPC_all<-read_csv("data-processed/TT_Schoolfield_TPC_all.csv")

GAM_ks_all<-read_csv("data-processed/GAM_ks_all.csv")
nitrate_decay_results<-read.csv(file="data-processed/nitrate_decay_rep_r-star.csv")
umax_predictions<-read_csv("data-processed/Norberg_model_umaxall_predictions.csv")


head(GAM_ks_all)
dim(umax_predictions)

#merge dataframes based on species, temperature, and runs
#prepare for merge
umax_predictions<-rename(umax_predictions, Temperature=temperature)
all_rstar<-merge(GAM_ks_all, umax_predictions, by=c("Species", "Temperature", "run"))
head(all_rstar)


# estimate R-star ------------
#R* = mKs / umax - m
#
#set m at 0.1
dim(umax_predictions)
dim(GAM_ks_all)
head(all_rstar)

fixed_ks<-1

all_rstar<-all_rstar %>%
  mutate(ksnozero=ifelse(ks==0, 0.06, ks)) %>%
  mutate(rstar=(ksnozero*0.1/ (all_rstar$umaxpredict-0.1))) %>% #calculate R-star for every simulation, m=0.1
  mutate(rstar_fixed_ks=(fixed_ks*0.1/ (all_rstar$umaxpredict-0.1))) %>%
  mutate(rstar_above_zero=ifelse(rstar<0, 0, rstar)) %>% #constrain r-star to 0
  mutate(rstar_fixed_ks_above_zero=ifelse(rstar_fixed_ks<0, 0, rstar_fixed_ks)) #constrain r-star to 0
         
#combine actual monod fits to predict rstar
all_umax<- all_ks_umax %>% filter(term=="umax") 
all_ks<- all_ks_umax %>% filter(term=="ks") 
rstar_estimate<- merge(all_ks, all_umax, by=c("Temperature", "run", "Species")) %>%
  rename(ks=estimate.x, umax=estimate.y) %>%
  mutate(rstar_calc=ks*0.1/(umax-0.1)) %>%
  mutate(rstar_calc_ks_fixed=fixed_ks*0.1/(umax-0.1))

all_rstar %>%
  ggplot(aes(y=umaxpredict, x=Temperature, group=run, color=as.factor(Species))) + facet_grid(~Species)+
  coord_cartesian(ylim = c(0, 10)) + geom_line(alpha=0.5) +
  theme_bw() + geom_point(data=rstar_estimate, aes(y=umax, x=Temperature), inherit.aes=FALSE, alpha=0.05)
ggsave("predicted_umax_with_temp_2015.pdf", width=7, height=2)

all_rstar %>%
  ggplot(aes(y=ks, x=Temperature, group=run, color=as.factor(Species))) + facet_grid(~Species)+
  coord_cartesian(ylim = c(0, 25)) + geom_line(alpha=0.5) +
  theme_bw() + geom_point(data=rstar_estimate, aes(y=ks, x=Temperature), inherit.aes=FALSE, alpha=0.05)
ggsave("predicted_ks_with_temp_2015.pdf", width=7, height=2)

all_rstar %>%
  ggplot(aes(y=rstar_above_zero, x=Temperature, group=run, color=as.factor(Species))) + facet_grid(~Species)+
  coord_cartesian(ylim = c(0, 5.5)) + geom_line(alpha=0.5)+
  theme_bw() + geom_point(data=rstar_estimate, aes(y=rstar_calc, x=Temperature), inherit.aes=FALSE, alpha=0.05)
ggsave("predicted_rstar_with_temp_2015.pdf", width=7, height=2)

all_rstar %>%
  ggplot(aes(y=rstar_fixed_ks_above_zero, x=Temperature, group=run, color=as.factor(Species))) + facet_grid(~Species)+
  coord_cartesian(ylim = c(0, 5.5)) + geom_line(alpha=0.5)+
  theme_bw() + geom_point(data=rstar_estimate, aes(y=rstar_calc_ks_fixed, x=Temperature), inherit.aes=FALSE, alpha=0.05)
ggsave("predicted_rstar_fixed_ks_with_temp_2015.pdf", width=7, height=2)

all_rstar_sum<-all_rstar %>%
  group_by(temperature, Species) %>%
  summarize(uci_rstar=quantile(rstar, 0.95), lci_rstar=quantile(rstar, 0.05), 
            median_rstar=quantile(rstar, 0.5))


all_rstar_sum %>%
  ggplot(aes(y=median_rstar, x=temperature)) + facet_grid(~Species)+
  geom_ribbon(aes(ymax=uci_rstar, ymin=lci_rstar, x=temperature), fill="grey70") +
  coord_cartesian(ylim = c(0, 10)) + geom_line() +
  geom_line(data=TT2017rstar, aes(y=K, x=temp), col="red")


