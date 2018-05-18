#goal: 
# fit the IDE model with Ks(g) [at end of model] simultaneously to all of the cell density
# this fulfils the "regression" approach
# and allows exploration of how Ks(g) changes with temperature 

library(broom)
library(purrr)
library(tidyverse)
library(minpack.lm)
library(nlstools)

TTfilteredN<-cbind(read_csv("data-processed/TTfilteredN.csv"), Species="TT")
CSfilteredN<-cbind(read_csv("data-processed/CSfilteredN.csv"), Species="CS")
CHfilteredN<-cbind(read_csv("data-processed/CHfilteredN.csv"), Species="CH")
ACfilteredN<-cbind(read_csv("data-processed/ACfilteredN.csv"), Species="AC")

#get fixed intercept for each species
TT_int<-TTfilteredN %>% 
  filter(DAY==0)  %>% 
  summarize(mean_init=mean(log.Particles.per.ml))
TTfilteredN<-TTfilteredN %>%
  mutate(mean_init=TT_int$mean_init)

CS_int<-CSfilteredN %>% 
  filter(DAY==0)  %>% 
  summarize(mean_init=mean(log.Particles.per.ml))
CSfilteredN<-CSfilteredN %>%
  mutate(mean_init=CS_int$mean_init)

CH_int<-CHfilteredN %>% 
  filter(DAY==0)  %>% 
  summarize(mean_init=mean(log.Particles.per.ml))
CHfilteredN<-CHfilteredN %>%
  mutate(mean_init=CH_int$mean_init)

AC_int<-ACfilteredN %>% 
  filter(DAY==0)  %>% 
  summarize(mean_init=mean(log.Particles.per.ml))
ACfilteredN<-ACfilteredN %>%
  mutate(mean_init=AC_int$mean_init)

#bring data together
alldata<-rbind(ACfilteredN, CHfilteredN, CSfilteredN, TTfilteredN)


#fitting IDE to cell density
#fixed k
all_species_IDE_kg<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks = 1), 
                lower=c(0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

write_csv(all_species_IDE_kg, "data-processed/all_species_IDE_kg")

#now take the root of this.
for (j in 1:length(unique(all_species_IDE_kg$Species))){ #for each species
  this_species<-unique(all_species_IDE_kg$Species)[j]
  b1<-filter(all_species_IDE_kg, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_kg, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_kg, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_kg, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_kg, Species==this_species)$estimate[5]
  ks<-filter(all_species_IDE_kg, Species==this_species)$estimate[6]
  
  roots_IDE<-data.frame()
  temp_IDE<-data.frame()
  
  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i]
    IDE_BR<-function(N){
      growth_rate<-((b1*exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature)))*(N/(ks + N)) -0.1
    }
    thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
    if(class(thisroot)!="try-error")
      roots_IDE<-rbind(roots_IDE, thisroot$root)
    if(class(thisroot)!="try-error")
      temp_IDE<-rbind(temp_IDE,Temperature)  
  }
  roots<-data.frame(root=roots_IDE[,1], temp=temp_IDE[,1], Species=this_species)
  assign(paste("roots",this_species,sep="_"), roots)
}


all_species_IDE_kg_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_IDE_kg_root, "data-processed/all_species_IDE_kg_root") #save

#plot the root
all_species_IDE_kg_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(all_species_IDE_kg$Species))){ #for each species
  this_species<-unique(all_species_IDE_kg$Species)[j]
  b1<-filter(all_species_IDE_kg, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_kg, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_kg, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_kg, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_kg, Species==this_species)$estimate[5]
  ks<-filter(all_species_IDE_kg, Species==this_species)$estimate[6]
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=ks, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

all_species_umax_kg<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_kg, "data-processed/all_species_umax_kg") #save

all_species_umax_kg<-read_csv("data-processed/all_species_umax_kg") #

#plot umax
all_species_umax_kg %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 1.5), xlim=c(0, 35)) +
  theme_bw() + ylab("Growth rate, r") 
ggsave("figures/growth_rate_2_spp.pdf", width=3, height=1.5)

#plot ks
all_species_umax_kg %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
IDE_all_species_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks = 1), 
                lower=c(0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                  ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                     -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks + N.Treatment)) * day,
                  data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks = 1), 
                  lower=c(0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 

#fitting IDE to cell density
#linear ks (=var)
all_species_IDE_kgvar<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01), 
                lower=c(0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

write_csv(all_species_IDE_kgvar, "data-processed/all_species_IDE_kgvar")


#now take the root of this.
for (j in 1:length(unique(all_species_IDE_kgvar$Species))){ #for each species
  this_species<-unique(all_species_IDE_kgvar$Species)[j]
  b1<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[5]
  ks_int<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[6]
  ks_slope<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[7]
  
  roots_IDE<-data.frame()
  temp_IDE<-data.frame()
  
  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i]
    IDE_BR<-function(N){
      growth_rate<-((b1*exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature)))*(N/(ks_int + ks_slope*Temperature + N)) -0.1
    }
    thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
    if(class(thisroot)!="try-error")
      roots_IDE<-rbind(roots_IDE, thisroot$root)
    if(class(thisroot)!="try-error")
      temp_IDE<-rbind(temp_IDE,Temperature)  
  }
  roots<-data.frame(root=roots_IDE[,1], temp=temp_IDE[,1], Species=this_species)
  assign(paste("roots",this_species,sep="_"), roots)
}


all_species_IDE_kg_var_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_IDE_kg_var_root, "data-processed/all_species_IDE_kg_var_root") #save

#plot the root
all_species_IDE_kg_var_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(all_species_IDE_kgvar$Species))){ #for each species
  this_species<-unique(all_species_IDE_kgvar$Species)[j]
  b1<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[5]
  ks_int<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[6]
  ks_slope<-filter(all_species_IDE_kgvar, Species==this_species)$estimate[7]
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=ks_int+ks_slope*newdata$Temperature, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

all_species_umax_kg_var<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_kg_var, "data-processed/all_species_umax_kg_var") #save

#plot umax
all_species_umax_kg_var %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
  theme_bw() + ylab("Predicted r-star, uM") 

#plot ks
all_species_umax_kg_var %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
IDE_all_species_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01), 
                lower=c(0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01), 
                lower=c(0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 

#
#
#
#
#fitting IDE to cell density
#polynomial ks (=poly)
all_species_IDE_kgpoly<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                lower=c(0,0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

write_csv(all_species_IDE_kgpoly, "data-processed/all_species_IDE_kgpoly")


#now take the root of this.
for (j in 1:length(unique(all_species_IDE_kgpoly$Species))){ #for each species
  this_species<-unique(all_species_IDE_kgpoly$Species)[j]
  b1<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[5]
  ks_int<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[6]
  ks_slope<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[7]
  so_ks<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[8]
  
  roots_IDE<-data.frame()
  temp_IDE<-data.frame()
  
  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i]
    IDE_BR<-function(N){
      growth_rate<-((b1*exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature)))*(N/(ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N)) -0.1
    }
    thisroot<-try(uniroot(IDE_BR, interval=c(0,400000000)), TRUE)
    if(class(thisroot)!="try-error")
      roots_IDE<-rbind(roots_IDE, thisroot$root)
    if(class(thisroot)!="try-error")
      temp_IDE<-rbind(temp_IDE,Temperature)  
  }
  roots<-data.frame(root=roots_IDE[,1], temp=temp_IDE[,1], Species=this_species)
  assign(paste("roots",this_species,sep="_"), roots)
}


all_species_IDE_kg_poly_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_IDE_kg_poly_root, "data-processed/all_species_IDE_kg_poly_root") #save

#plot the root
all_species_IDE_kg_poly_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(all_species_IDE_kgpoly$Species))){ #for each species
  this_species<-unique(all_species_IDE_kgpoly$Species)[j]
  b1<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[5]
  ks_int<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[6]
  ks_slope<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[7]
  so_ks<-filter(all_species_IDE_kgpoly, Species==this_species)$estimate[8]
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=ks_int+ks_slope*newdata$Temperature+so_ks*newdata$Temperature^2, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

all_species_umax_kg_poly<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_kg_poly, "data-processed/all_species_umax_kg_poly") #save

#plot umax
all_species_umax_kg_poly %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
  theme_bw() + ylab("Predicted r-star, uM") 

#plot ks
all_species_umax_kg_poly %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
IDE_all_species_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                lower=c(0,0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N.Treatment)) * day,
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                lower=c(0,0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 



#plot all the umax's
all_species_umax_kg_poly %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
  geom_line(data=all_species_umax_kg_var, linetype=2) +
  geom_line(data=all_species_umax_kg, linetype=3) +
  theme_bw() + ylab("Predicted r-star, uM") 


#plot all the roots
all_species_IDE_kg_poly_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 100)) + 
  geom_line(data=all_species_IDE_kg_var_root, linetype=2) +
  geom_line(data=all_species_IDE_kg_root, linetype=3) +
  theme_bw() + ylab("Predicted r-star, uM") 

ggsave("r_star_predictions_IDE_kg.pdf", width=7, height=2)
