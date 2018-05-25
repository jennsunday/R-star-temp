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


#bring data together
alldata<-rbind(ACfilteredN, CHfilteredN, CSfilteredN, TTfilteredN) 


#fitting IDE to cell density
#fixed k
#nlsmultstart
IDE_kb<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (ks + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.024, b2= 0.14, d0=0, d1=0, d2=3,  ks = 0.1), 
                lower=c(0,0,-5,0,0,-5,0), control=nls.control(maxiter = 1000))))
#View(IDE_kb)
write_csv(IDE_kb, "data-processed/IDE_kb")


#now take the root of this.
newdata<-data.frame(Temperature=seq(3, 42, 0.1))

for (j in 1:length(unique(IDE_kb$Species))){ #for each species
  this_species<-unique(IDE_kb$Species)[j]
  b1<-filter(IDE_kb, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kb, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kb, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kb, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kb, Species==this_species & term=="d2")$estimate
  ks<-filter(IDE_kb, Species==this_species & term=="ks")$estimate

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


IDE_kb_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(IDE_kb_root, "data-processed/IDE_kb_root") #save

#plot the root
IDE_kb_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 10)) + 
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(IDE_kb$Species))){ #for each species
  this_species<-unique(IDE_kb$Species)[j]
  b1<-filter(IDE_kb, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kb, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kb, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kb, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kb, Species==this_species & term=="d2")$estimate
  ks<-filter(IDE_kb, Species==this_species & term=="ks")$estimate
  
  
umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                    ks_temp=ks, Species=this_species, temp=newdata$Temperature)
assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

IDE_kb_umax_ks<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(IDE_kb_umax_ks, "data-processed/IDE_kb_umax_ks") #save

#plot umax
IDE_kb_umax_ks %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
  theme_bw() + ylab("umax") 

#plot ks
IDE_kb_umax_ks %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
  theme_bw() + ylab("ks_birth") 

#glance model summary
IDE_kb_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (ks + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.024, b2= 0.14, d0=0, d1=0, d2=3,  ks = 0.1), 
                lower=c(0,0,-5,0,0,-5,0), control=nls.control(maxiter = 1000))))

#plot predicted by fitted
alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                  ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (ks + N.Treatment)))
                                     -(d0 + d1*exp(d2*Temperature)))) * day,
                  data= .,  start=list(N_init=5.7, b1 = 0.024, b2= 0.14, d0=0, d1=0, d2=3,  ks = 0.1), 
                  lower=c(0,0,-5,-0.05,0,-5,0), control=nls.control(maxiter = 1000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
theme_bw() + ylab("predicted") 

#
#letting ks change with temp - linear
#
IDE_kb_var<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(N_init=5.7, b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05, int_ks = 0, slope_ks=0.1), 
                lower=c(0,0,-5,0,0,-5,0,0), control=nls.control(maxiter = 1000))))

write_csv(IDE_kb_var, "data-processed/IDE_kb_var")
#View(IDE_kb_var)
#
#
#

#now take the root of this.
newdata<-data.frame(Temperature=seq(3, 38, 0.1))

for (j in 1:length(unique(IDE_kb_var$Species))){ #for each species
  this_species<-unique(IDE_kb_var$Species)[j]
  b1<-filter(IDE_kb_var, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kb_var, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kb_var, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kb_var, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kb_var, Species==this_species & term=="d2")$estimate
  int_ks<-filter(IDE_kb_var, Species==this_species & term=="int_ks")$estimate
  slope_ks<-filter(IDE_kb_var, Species==this_species & term=="slope_ks")$estimate

roots_IDE<-data.frame()
temp_IDE<-data.frame()

for(i in 1:length(newdata$Temperature)){
  Temperature<-newdata$Temperature[i]
  IDE_BR<-function(N){
    growth_rate<-b1*exp(b2*Temperature)*(N/(int_ks+(slope_ks*Temperature)+ N))-
      (d0 + d1*exp(d2*Temperature))-0.1
  }
  thisroot<-try(uniroot(IDE_BR, interval=c(0,4000000000)), TRUE)
  if(class(thisroot)!="try-error")
    roots_IDE<-rbind(roots_IDE, thisroot$root)
  if(class(thisroot)!="try-error")
    temp_IDE<-rbind(temp_IDE,Temperature)  
}
roots<-data.frame(root=roots_IDE[,1], temp=temp_IDE[,1], Species=this_species)
assign(paste("roots",this_species,sep="_"), roots)
}

IDE_kb_var_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(IDE_kb_var_root, "data-processed/IDE_kb_var_root") #save


IDE_kb_var_root<-read_csv("data-processed/IDE_kb_var_root")

IDE_kb_var_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 10)) + theme_bw()

#plot ks and umax parameters
for (j in 1:length(unique(IDE_kb_var$Species))){ #for each species
  this_species<-unique(IDE_kb_var$Species)[j]
  b1<-filter(IDE_kb_var, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kb_var, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kb_var, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kb_var, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kb_var, Species==this_species & term=="d2")$estimate
  int_ks<-filter(IDE_kb_var, Species==this_species & term=="int_ks")$estimate
  slope_ks<-filter(IDE_kb_var, Species==this_species & term=="slope_ks")$estimate
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=int_ks+slope_ks*newdata$Temperature, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

umax_kb_var<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(umax_kb_var, "data-processed/umax_kb_var") #save

#plot umax
umax_kb_var %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
theme_bw() + ylab("Predicted umax, uM") 


#plot ks
umax_kb_var %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
theme_bw() + ylab("Ks") 

#glance model summary
IDE_kb_var_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(N_init=5.7, b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05, int_ks = -2, slope_ks=0.1), 
                lower=c(0,0,-5,0,0,-5,0,0), control=nls.control(maxiter = 1000))))

#plot predicted by fitted
alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                  ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + N.Treatment)))
                                     -(d0 + d1*exp(d2*Temperature)))) * day,
                  data= .,  start=list(N_init=5.7, b1 = 11, b2= 0.1, d0=11, d1=0.14, d2=0.05, int_ks = 0, slope_ks=0.1), 
                  lower=c(0,0,-5,0,0,-5,0,0), control=nls.control(maxiter = 1000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("modeled") 


#
#letting ks change with temp - polynomial
#
IDE_kb_poly<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + int_ks+so_ks*Temperature^2+ N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.024, b2= 0.14, d0=0, d1=0, d2=3,  int_ks = 0, slope_ks=0.01, so_ks=0.001), 
                lower=c(0,0,0,-4,0,0,-4,-20,0), control=nls.control(maxiter = 1000)))) 

write_csv(IDE_kb_poly, "data-processed/IDE_kb_poly")
View(IDE_kb_poly)
#
#
#

#now take the root of this.
for (j in 1:length(unique(IDE_kb_poly$Species))){ #for each species
  this_species<-unique(IDE_kb_poly$Species)[j]
  b1<-filter(IDE_kb_poly, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kb_poly, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kb_poly, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kb_poly, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kb_poly, Species==this_species & term=="d2")$estimate
  int_ks<-filter(IDE_kb_poly, Species==this_species & term=="int_ks")$estimate
  slope_ks<-filter(IDE_kb_poly, Species==this_species & term=="slope_ks")$estimate
  so_ks<-filter(IDE_kb_poly, Species==this_species & term=="so_ks")$estimate
  
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

IDE_kb_poly_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(IDE_kb_poly_root, "data-processed/IDE_kb_poly_root") #save


#plot the root
IDE_kb_poly_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 100)) + 
  theme_bw() + ylab("Predicted r-star, uM") 

#get ks and umax
for (j in 1:length(unique(IDE_kb_poly$Species))){ #for each species
  this_species<-unique(IDE_kb_poly$Species)[j]
  b1<-filter(IDE_kb_poly, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kb_poly, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kb_poly, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kb_poly, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kb_poly, Species==this_species & term=="d2")$estimate
  int_ks<-filter(IDE_kb_poly, Species==this_species & term=="int_ks")$estimate
  slope_ks<-filter(IDE_kb_poly, Species==this_species & term=="slope_ks")$estimate
  so_ks<-filter(IDE_kb_poly, Species==this_species & term=="so_ks")$estimate

  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=int_ks+int_slope*newdata$Temperature+so_ks*(newdata$Temperature)^2, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

umax_kb_poly<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(umax_kb_poly, "data-processed/umax_kb_poly") #save

#plot umax
umax_kb_poly %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species)
theme_bw() + ylab("Predicted r-star, uM") 
ggsave("figures/umax_kb_poly", width=8, height=4)

#plot ks
umax_kb_poly %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species)
theme_bw() + ylab("Predicted r-star, uM") 


#now check fit
IDE_kb_poly_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + int_ks+so_ks*Temperature^2+ N.Treatment)))
                                   -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.024, b2= 0.14, d0=0, d1=0, d2=3, int_ks = 0, slope_ks=0.01, so_ks=0.001), 
                lower=c(0,0,0,-4,0,0,-4,-20,0), control=nls.control(maxiter = 1000)))) 


alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature) * (N.Treatment / (int_ks+slope_ks*Temperature + int_ks+so_ks*Temperature^2+ N.Treatment)))
                             -(d0 + d1*exp(d2*Temperature)))) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.024, b2= 0.14, d0=0, d1=0, d2=3, int_ks = 0, slope_ks=0.01, so_ks=0.001), 
                lower=c(0,0,0,-4,0,0,-4,-20,0), control=nls.control(maxiter = 1000))))  %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("predicted") 


