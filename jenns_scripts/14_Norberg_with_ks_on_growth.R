### fit TPC to 2015 umax data 
#goal:fit Norberg model to umax data for 1000 bootstraps for each of 4 species
#worried this will take a long time
#so reduce the "searching" scope so that the time to fit each model is shorter.

library(broom)
library(purrr)
library(tidyverse)
library(minpack.lm)
library(nlstools)

#read in data
TTfilteredN<-cbind(read_csv("data-processed/TTfilteredN.csv"), Species="TT")
CSfilteredN<-cbind(read_csv("data-processed/CSfilteredN.csv"), Species="CS")
CHfilteredN<-cbind(read_csv("data-processed/CHfilteredN.csv"), Species="CH")
ACfilteredN<-cbind(read_csv("data-processed/ACfilteredN.csv"), Species="AC")


#retrieve potential starting values from 
Norberg_TPC_fits_2015

####################################
#fitting NB to cell density
#fixed k
####################################
all_species_NB_kg<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks + N.Treatment)) * day,
                data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks = 1), 
                lower=c(0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))
write_csv(all_species_NB_kg, "data-processed/all_species_NB_kg")

#now take the root of this.
for (j in 1:length(unique(all_species_NB_kg$Species))){ #for each species
  this_species<-unique(all_species_NB_kg$Species)[j]
  a<-filter(all_species_NB_kg, Species==this_species)$estimate[1]
  b<-filter(all_species_NB_kg, Species==this_species)$estimate[2]
  z<-filter(all_species_NB_kg, Species==this_species)$estimate[3]
  w<-filter(all_species_NB_kg, Species==this_species)$estimate[4]
  ks<-filter(all_species_NB_kg, Species==this_species)$estimate[5]

  roots_NB<-data.frame() #make an empty dataframe
  temp_NB<-data.frame() #make an empty dataframe

  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i] #for every half-degree from 3 to 38
    NB_BR<-function(N){
      growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2) * N/(ks + N) - 0.1}
    thisroot<-try(uniroot(NB_BR, interval=c(0,400000000)), TRUE) #try taking the root
    if(class(thisroot)!="try-error") #if you don't get an error,
      roots_NB<-rbind(roots_NB, thisroot$root) #add this to the output
    if(class(thisroot)!="try-error") #if you don't get an error
      temp_NB<-rbind(temp_NB,Temperature) #mark which temperature
  }
  
  roots<-data.frame(root=roots_NB[,1], temp=temp_NB[,1], Species=this_species)
  assign(paste("roots",this_species,sep="_"), roots)
}

all_species_NB_kg_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_NB_kg_root, "data-processed/all_species_NB_kg_root") #save

#plot the root
all_species_NB_kg_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(all_species_NB_kg$Species))){ #for each species
  this_species<-unique(all_species_NB_kg$Species)[j]
  a<-filter(all_species_NB_kg, Species==this_species)$estimate[1]
  b<-filter(all_species_NB_kg, Species==this_species)$estimate[2]
  z<-filter(all_species_NB_kg, Species==this_species)$estimate[3]
  w<-filter(all_species_NB_kg, Species==this_species)$estimate[4]
  ks<-filter(all_species_NB_kg, Species==this_species)$estimate[5]
  
  umax_ks<-data.frame(umax_temp=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                      ks_temp=ks, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

all_species_umax_NB_kg<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_NB_kg, "data-processed/all_species_umax_NB_kg") #save

all_species_umax_NB_kg<-read_csv("data-processed/all_species_umax_NB_kg") #

#plot umax
all_species_umax_NB_kg %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 1.5), xlim=c(0, 35)) +
  theme_bw() + ylab("Growth rate, r") 
ggsave("figures/growth_rate_NB_kg.pdf", width=3, height=1.5)

#plot ks
all_species_umax_NB_kg %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
NB_kg_all_species_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                  ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks + N.Treatment)) * day,
                  data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks = 1), 
                  lower=c(0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                   ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks + N.Treatment)) * day,
                   data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks = 1), 
                   lower=c(0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 


####################################
#fitting NB to cell density
#linear relationship between k and temperature
####################################

all_species_NB_kg_var<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks_int = 1, ks_slope=0.01), 
                lower=c(0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))
write_csv(all_species_NB_kg_var, "data-processed/all_species_NB_kg_var")

#now take the root of this.
for (j in 1:length(unique(all_species_NB_kg_var$Species))){ #for each species
  this_species<-unique(all_species_NB_kg_var$Species)[j]
  a<-filter(all_species_NB_kg_var, Species==this_species)$estimate[1]
  b<-filter(all_species_NB_kg_var, Species==this_species)$estimate[2]
  z<-filter(all_species_NB_kg_var, Species==this_species)$estimate[3]
  w<-filter(all_species_NB_kg_var, Species==this_species)$estimate[4]
  ks_int<-filter(all_species_NB_kg_var, Species==this_species)$estimate[5]
  ks_slope<-filter(all_species_NB_kg_var, Species==this_species)$estimate[6]
  
  roots_NB_ks_var<-data.frame() #make an empty dataframe
  temp_NB_ks_var<-data.frame() #make an empty dataframe
  
  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i] #for every half-degree from 3 to 38
    NB_BR<-function(N){
      growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2) * N/(int_ks+slope_ks*Temperature + N) - 0.1}
    thisroot<-try(uniroot(NB_BR, interval=c(0,400000000)), TRUE) #try taking the root
    if(class(thisroot)!="try-error") #if you don't get an error,
      roots_NB_ks_var<-rbind(roots_NB_ks_var, thisroot$root) #add this to the output
    if(class(thisroot)!="try-error") #if you don't get an error
      temp_NB_ks_var<-rbind(temp_NB_ks_var,Temperature) #mark which temperature
  }
  
  roots<-data.frame(root=roots_NB_ks_var[,1], temp=temp_NB_ks_var[,1], Species=this_species)
  assign(paste("roots",this_species,sep="_"), roots)
}

all_species_NB_kg_var_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_NB_kg_var_root, "data-processed/all_species_NB_kg_var_root") #save

#plot the root
all_species_NB_kg_var_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(all_species_NB_kg_var$Species))){ #for each species
  this_species<-unique(all_species_NB_kg_var$Species)[j] #get the name of the species
  a<-filter(all_species_NB_kg_var, Species==this_species)$estimate[1]
  b<-filter(all_species_NB_kg_var, Species==this_species)$estimate[2]
  z<-filter(all_species_NB_kg_var, Species==this_species)$estimate[3]
  w<-filter(all_species_NB_kg_var, Species==this_species)$estimate[4]
  ks_int<-filter(all_species_NB_kg_var, Species==this_species)$estimate[5]
  ks_slope<-filter(all_species_NB_kg_var, Species==this_species)$estimate[6]
  
  umax_ks<-data.frame(umax_temp=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                      ks_temp=ks_int+ks_slope*newdata$Temperature, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}


all_species_umax_NB_kg_var<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_NB_kg_var, "data-processed/all_species_umax_NB_kg_var") #save

all_species_umax_NB_kg_var<-read_csv("data-processed/all_species_umax_NB_kg_var") #

#plot umax
all_species_umax_NB_kg_var %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 1.5), xlim=c(0, 35)) +
  theme_bw() + ylab("Growth rate, r") 
ggsave("figures/growth_rate_NB_kg_var.pdf", width=3, height=1.5)

#plot ks
all_species_umax_NB_kg_var %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
NB_kg_var_all_species_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                  ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                  data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks_int = 1, ks_slope=0.01), 
                  lower=c(0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                   ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                   data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks_int = 1, ks_slope=0.01), 
                   lower=c(0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 


####################################
#fitting NB to cell density
#polynomial relationship between k and temperature
####################################

all_species_NB_kg_poly<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2  + N.Treatment)) * day,
                data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                lower=c(0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))
write_csv(all_species_NB_kg_poly, "data-processed/all_species_NB_kg_poly")

#now take the root of this.
for (j in 1:length(unique(all_species_NB_kg_poly$Species))){ #for each species
  this_species<-unique(all_species_NB_kg_poly$Species)[j]
  a<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[1]
  b<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[2]
  z<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[3]
  w<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[4]
  ks_int<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[5]
  ks_slope<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[6]
  so_ks<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[7]
  
  roots_NB_ks_poly<-data.frame() #make an empty dataframe
  temp_NB_ks_poly<-data.frame() #make an empty dataframe
  
  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i] #for every half-degree from 3 to 38
    NB_BR<-function(N){
      growth_rate<-a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2) * N/(ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N) - 0.1}
    thisroot<-try(uniroot(NB_BR, interval=c(0,400000000)), TRUE) #try taking the root
    if(class(thisroot)!="try-error") #if you don't get an error,
      roots_NB_ks_poly<-rbind(roots_NB_ks_poly, thisroot$root) #add this to the output
    if(class(thisroot)!="try-error") #if you don't get an error
      temp_NB_ks_poly<-rbind(temp_NB_ks_poly,Temperature) #mark which temperature
  }
  
  roots<-data.frame(root=roots_NB_ks_poly[,1], temp=temp_NB_ks_poly[,1], Species=this_species)
  assign(paste("roots",this_species,sep="_"), roots)
}

all_species_NB_kg_poly_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_NB_kg_poly_root, "data-processed/all_species_NB_kg_poly_root") #save

#plot the root
all_species_NB_kg_poly_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(all_species_NB_kg_poly$Species))){ #for each species
  this_species<-unique(all_species_NB_kg_poly$Species)[j] #get the name of the species
  a<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[1]
  b<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[2]
  z<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[3]
  w<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[4]
  ks_int<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[5]
  ks_slope<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[6]
  so_ks<-filter(all_species_NB_kg_poly, Species==this_species)$estimate[7]
  
  umax_ks<-data.frame(umax_temp=a*exp(b*newdata$Temperature)*(1-((newdata$Temperature-z)/(w/2))^2), 
                      ks_temp=ks_int+ks_slope*newdata$Temperature+so_ks*newdata$Temperature^2, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}


all_species_umax_NB_kg_poly<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_NB_kg_poly, "data-processed/all_species_umax_NB_kg_poly") #save

all_species_umax_NB_kg_poly<-read_csv("data-processed/all_species_umax_NB_kg_poly") #

#plot umax
all_species_umax_NB_kg_poly %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 1.5), xlim=c(0, 35)) +
  theme_bw() + ylab("Growth rate, r") 
ggsave("figures/growth_rate_NB_kg_var.pdf", width=3, height=1.5)

#plot ks
all_species_umax_NB_kg_poly %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
NB_kg_poly_all_species_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                  ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2  + N.Treatment)) * day,
                  data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                  lower=c(0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                   ~ mean_init[1] + (a*exp(b*Temperature)*(1-((Temperature-z)/(w/2))^2)) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2  + N.Treatment)) * day,
                   data= .,  start=list(a = 0.30661082, b=0.05808612, z=20.310270, w=26.58936, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                   lower=c(0,0,0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 
