#goal: 
# fit the whole IDE model with Ks(b) simultaneously to all of the cell density
# this fulfils the "regression" approach
# and allows exploration of how Ks(b) changes with temperature 

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
all_species_IDE<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (ks + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(b1 = 1, b2= 0.1, d0=1, d1=0.14, d2=0.05,  ks = 0.1), 
                lower=c(0,0,0,0,0,0))))

write_csv(all_species_IDE, "data-processed/all_species_IDE")

#now take the root of this.
for (j in 1:length(unique(all_species_IDE$Species))){ #for each species
  this_species<-unique(all_species_IDE$Species)[j]
  all_species_IDE
  b1<-filter(all_species_IDE, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE, Species==this_species)$estimate[5]
  ks<-filter(all_species_IDE, Species==this_species)$estimate[6]

  roots_IDE<-data.frame()
  temp_IDE<-data.frame()

  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i]
    IDE_BR<-function(N){
      growth_rate<-b1*exp(b2*Temperature)*(N/(ks + N))-(d0 + d1*exp(d2*Temperature)) -0.1
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


all_species_IDE_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_IDE_root, "data-processed/all_species_IDE_root") #save

#plot the root
all_species_IDE_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 20)) + 
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(all_species_IDE$Species))){ #for each species
  this_species<-unique(all_species_IDE$Species)[j]
  all_species_IDE
  b1<-filter(all_species_IDE, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE, Species==this_species)$estimate[5]
  ks<-filter(all_species_IDE, Species==this_species)$estimate[6]
  
umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                    ks_temp=ks, Species=this_species, temp=newdata$Temperature)
assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

all_species_umax_ks<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_ks, "data-processed/all_species_umax_ks") #save

#plot umax
all_species_umax_ks %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
  theme_bw() + ylab("Predicted r-star, uM") 

#plot ks
all_species_umax_ks %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
  theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
IDE_all_species_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (ks + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  ks = 0.1), 
                lower=c(0,0,0,0,0,0))))

#plot observed by fitted
alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                  ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (ks + N.Treatment)))
                                     -(d0 + d1*exp(d2*Temperature)))) * day,
                  data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  ks = 0.1), 
                  lower=c(0,0,0,0,0,0)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
theme_bw() + ylab("observed") 

#
#letting ks change with temp - linear
#
all_species_IDE_ks_var<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  int_ks = -2, slope_ks=0.1), 
                lower=c(0,0,0,0,0,0,0))))

write_csv(all_species_IDE_ks_var, "data-processed/all_species_IDE_ks_var")

#
#
#

#now take the root of this.
for (j in 1:length(unique(all_species_IDE_ks_var$Species))){ #for each species
  this_species<-unique(all_species_IDE_ks_var$Species)[j]
  b1<-filter(all_species_IDE_var, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_var, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_var, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_var, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_var, Species==this_species)$estimate[5]
  int_ks<-filter(all_species_IDE_var, Species==this_species)$estimate[6]
  slope_ks<-filter(all_species_IDE_var, Species==this_species)$estimate[7]

roots_IDE<-data.frame()
temp_IDE<-data.frame()

for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i]
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(int_ks+slope_ks*Temperature + N))-
      (d0 + d1*exp(d2*Temperature))-0.1
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

all_species_IDE_root_ks_var<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_IDE_root, "data-processed/all_species_IDE_root_ks_var") #save


all_species_IDE_root_var<-read_csv("data-processed/all_species_IDE_root_var")
all_species_IDE_root<-read_csv("data-processed/all_species_IDE_root")

all_species_IDE_root_ks_var %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 100)) + theme_bw() +
  geom_line(data=all_species_IDE_root, linetype=2)

#plot ks and umax parameters
for (j in 1:length(unique(all_species_IDE_ks_var$Species))){ #for each species
  this_species<-unique(all_species_IDE_ks_var$Species)[j]
  b1<-filter(all_species_IDE_ks_var, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_ks_var, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_ks_var, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_ks_var, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_ks_var, Species==this_species)$estimate[5]
  int_ks <-filter(all_species_IDE_ks_var, Species==this_species)$estimate[6]
  int_slope <-filter(all_species_IDE_ks_var, Species==this_species)$estimate[7]
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=int_ks+int_slope*newdata$Temperature, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

all_species_umax_ks_var<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_ks_var, "data-processed/all_species_umax_ks_var") #save

#plot umax
all_species_umax_ks_var %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
theme_bw() + ylab("Predicted umax, uM") 


#plot ks
all_species_umax_ks_var %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
all_species_IDE_ks_var_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  int_ks = -2, slope_ks=0.1), 
                lower=c(0,0,0,0,0,0,0))))

#plot observed by fitted
alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                  ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + N.Treatment)))
                                     -(d0 + d1*exp(d2*Temperature)))) * day,
                  data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  int_ks = -2, slope_ks=0.1), 
                  lower=c(0,0,0,0,0,0,0)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 



#
#letting ks change with temp - polynomial
#
all_species_IDE_ks_poly<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + int_ks+so_ks*Temperature^2+ N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  int_ks = 0, slope_ks=0.01, so_ks=0.001), 
                lower=c(0,0,0,0,0,0,0,0))))

write_csv(all_species_IDE_ks_poly, "data-processed/all_species_IDE_ks_poly")

#
#
#

#now take the root of this.
for (j in 1:length(unique(all_species_IDE_ks_poly$Species))){ #for each species
  this_species<-unique(all_species_IDE_ks_poly$Species)[j]
  b1<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[5]
  int_ks<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[6]
  slope_ks<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[7]
  so_ks<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[8]
  
  roots_IDE<-data.frame()
  temp_IDE<-data.frame()
  
  for(i in 1:length(newdata$Temperature)){
    Temperature<-newdata$Temperature[i]
    IDE_BR<-function(N){
      growth_rate<-b1*exp(b2*Temperature)*(N/(int_ks+slope_ks*Temperature + int_ks+so_ks*Temperature^2+ N))-
        (d0 + d1*exp(d2*Temperature))-0.1
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

all_species_IDE_root_ks_poly<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(all_species_IDE_root_ks_poly, "data-processed/all_species_IDE_root_ks_poly") #save


#plot the root
all_species_IDE_root_ks_poly %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 100)) + 
  theme_bw() + ylab("Predicted r-star, uM") 

#get ks and umax
for (j in 1:length(unique(all_species_IDE_ks_poly$Species))){ #for each species
  this_species<-unique(all_species_IDE_ks_poly$Species)[j]
  b1<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[1]
  b2<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[2]
  d0<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[3]
  d1<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[4]
  d2<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[5]
  int_ks <-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[6]
  int_slope <-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[7]
  so_ks<-filter(all_species_IDE_ks_poly, Species==this_species)$estimate[8]

  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=int_ks+int_slope*newdata$Temperature+so_ks*(newdata$Temperature)^2, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

all_species_umax_ks_poly<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(all_species_umax_ks_poly, "data-processed/all_species_umax_ks_poly") #save

#plot umax
all_species_umax_ks_poly %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species)
theme_bw() + ylab("Predicted r-star, uM") 
ggsave("figures/umax_ks_poly.pdf", width=8, height=4)

#plot ks
all_species_umax_ks_poly %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species)
theme_bw() + ylab("Predicted r-star, uM") 


#now check fit
all_species_IDE_ks_poly_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + int_ks+so_ks*Temperature^2+ N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  int_ks = 0, slope_ks=0.01, so_ks=0.001), 
                lower=c(0,0,0,0,0,0,0,0))))


all_species_IDE_ks_poly_pred<-alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                ~ mean_init[1] + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + int_ks+so_ks*Temperature^2+ N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05,  int_ks = 0, slope_ks=0.01, so_ks=0.001), 
                lower=c(0,0,0,0,0,0,0,0))))


all_species_IDE_ks_poly_pred %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 

#Compare AIC's of each model
IDE_all_species_glance$AIC # fixed Ks is best for species TT, species CS
all_species_IDE_ks_var_glance$AIC #linear Ks is best for species AC
all_species_IDE_ks_poly_glance$AIC #polynomial Ks best for species CH


#plot all the roots together
all_species_IDE_root_ks_poly<-read_csv("data-processed/all_species_IDE_root_ks_poly")
all_species_IDE_root_var<-read_csv("data-processed/all_species_IDE_root_var")
all_species_IDE_root<-read_csv("data-processed/all_species_IDE_root")

all_species_IDE_root_ks_poly %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 100)) + 
  geom_line(data=all_species_IDE_root_var, linetype=2) +
  geom_line(data=all_species_IDE_root, linetype=3) +
  theme_bw() + ylab("Predicted r-star, uM") 