#get root for every temperature
#establish temperature frame

library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(cowplot)

newdata<-data.frame(Temperature=seq(-3, 50, 0.5))


# Species 1
a<-0.102
b<-0.143
z<-(-20)
w<-107.3
ks<-3

roots_NB<-data.frame() #make an empty dataframe - temporary
temp_NB<-data.frame() #make an empty dataframe - temporary

for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i] #for every half-degree from 1 to 50
  NB_fix<-function(N){
    growth_rate<-(a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2) * (N)/(ks + N)) - 0.1} # don't include tn so that actual root can be calculated
  thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
  roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                       Temperature=Temperature)) #add this to the output, setting root at a large number (200) if not within interval
}
roots_ks_umax<-mutate(roots_NB, umax=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                      ks=ks, Species=1)

# Species 2
a<-0.04
b<-0.14
z<-(-10)
w<-107
ks<-1
roots_NB<-data.frame() #make an empty dataframe - temporary
temp_NB<-data.frame() #make an empty dataframe - temporary

for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i] #for every half-degree from 1 to 50
  NB_fix<-function(N){
    growth_rate<-(a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2) * (N)/(ks + N)) - 0.1} # don't include tn so that actual root can be calculated
  thisroot<-try(uniroot(NB_fix, interval=c(0,200)), TRUE) #try taking the root
  roots_NB<-rbind(roots_NB, data.frame(root=ifelse(class(thisroot)!="try-error", thisroot$root, 200),
                                       Temperature=Temperature)) #add this to the output, setting root at a large number (200) if not within interval
}
roots_ks_umax2<-mutate(roots_NB, umax=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                      ks=ks, Species=2)

all_roots_ks_umax<-rbind(roots_ks_umax, roots_ks_umax2) 

all_roots_ks_umax<-all_roots_ks_umax %>%
  group_by(Temperature) %>%
  mutate(winner=Species[which.min(root)]) %>%
  ungroup()

View(all_roots_ks_umax)

P1<-all_roots_ks_umax %>%
  ggplot(aes(colour=as.factor(Species))) + geom_line(aes(x=Temperature, y=umax)) +
  coord_cartesian(ylim = c(0, 2), xlim = c(-3, 50))

P2<-all_roots_ks_umax %>%
  ggplot(aes(colour=as.factor(Species))) + geom_line(aes(x=Temperature, y=root)) +
  coord_cartesian(ylim = c(0, 5), xlim = c(-3, 50))

View(all_roots_ks_umax %>%
  filter(winner==1))

P3<-all_roots_ks_umax %>%
  filter(root<199) %>%
  ggplot() + geom_line(aes(x=Temperature, y=1, colour=as.factor(winner)), size=1) +
  coord_cartesian(ylim = c(0, 5), xlim = c(-3, 50))

together<-grid.arrange(P1, P2, P3,  nrow = 3)
ggsave("figures/didactic.pdf", together, width=5, height=7)

