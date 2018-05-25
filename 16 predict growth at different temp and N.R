all_roots_ks_umax_boot<-read_csv("data-processed/all_roots_ks_umax_boot_NB_var.csv")
all_roots_ks_umax<-read_csv("data-processed/all_roots_ks_umax_NB_var.csv")

head(all_roots_ks_umax)

### fit TPC to 2015 umax data 
#goal:fit Norberg model to umax data for 20 bootstraps for each of 4 species

newdata<-data.frame(Temperature=seq(-3, 50, 0.5))

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
  
  #get growth rate for 5 Ns
  test<-function (N,Temp) {
    a*exp(b*Temp)*(1-((Temp-z)/(w/2))^2) * N/(ks + N)
    }

df<-data.frame(N=seq(0.5, 40, length.out=25),
               Temp=seq(-3, 50, , length.out=25))

predict_growth<-df %>%
  expand(N, Temp) %>%
  mutate(growth=test(N,Temp))

predict_growth %>%
  ggplot(aes(x=Temp, y=growth, colour=as.factor(N))) + geom_line() +
  coord_cartesian(ylim = c(0, 1.5), xlim = c(0, 50))
  
test(test_vars)

  sum.of.squares <- function(x,y) {
    x^2 + y^2
  }
  sum.of.squares(3,4)
  
    thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
    roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                         Temperature=Temperature)) #add this to the output, setting root at a large number (200) if not within interval
  }
  roots_ks_umax<-mutate(roots_NB, umax=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                        ks=ks, Species=unique(alldata$Species)[k])
  all_roots_ks_umax<-rbind(all_roots_ks_umax, roots_ks_umax)
}
