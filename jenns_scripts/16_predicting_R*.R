

#once I settle on a plan for curve-fitting, I need to combine both curve fits within the same 
#draw, because ks and umax are *correlated* within runs. But for now just curve fitting sorted.
#combining draws will probably reduce extremes.

TT_Schoolfield_TPC_all<-read_csv("data-processed/TT_Schoolfield_TPC_all.csv")
GAM_ks_all<-read_csv("data-processed/GAM_ks_all.csv")
obs_rstar_2017<-read_csv("data-processed/logistic_N_decay_fits_r-star.csv")
TT2017rstar<-obs_rstar_2017 %>% filter(species==1)

plot(GAM_ks_all$ks~GAM_ks_all$Temperature, type="l")
plot(TT_Schoolfield_TPC_all$umax~TT_Schoolfield_TPC_all$Temperature, type="l")

# estimate R-star ------------
#R* = mKs / umax - m
#
#set m at 0.1
dim(TT_Schoolfield_TPC_all)
dim(GAM_ks_all)

all_rstar<-cbind(TT_Schoolfield_TPC_all, GAM_ks_all)
all_rstar<-all_rstar[,1:5] #get rid of duplicated columns
all_rstar$rstar<-10*0.1/ (all_rstar$umax-0.1)#calculate R-star for every simulation, m=0.1
all_rstar$rstar<-ifelse(all_rstar$rstar<0, 0, all_rstar$rstar)#constrain r-star to 0


all_rstar_sum<-all_rstar %>%
  group_by(Temperature) %>%
  summarize(uci_rstar=quantile(rstar, 0.95), lci_rstar=quantile(rstar, 0.05), 
            median_rstar=quantile(rstar, 0.5))

head(all_rstar_sum)
unique(all_rstar$run)

plot(all_rstar_sum$median_rstar~all_rstar_sum$Temperature, type="l", ylim=c(0,10))
points(all_rstar$rstar~all_rstar$Temperature, col="#00000030")

unique(all_rstar_sum$Temperature)
all_rstar_sum %>%
  ggplot(aes(y=median_rstar, x=Temperature)) + 
  geom_ribbon(aes(ymax=uci_rstar, ymin=lci_rstar, x=Temperature), fill="grey70") +
  coord_cartesian(ylim = c(0, 10)) + geom_line() +
  geom_line(data=TT2017rstar, aes(y=K, x=temp), col="red")


