#goal: 
# fit the IDE model with Ks(g) [at end of model] simultaneously to all of the cell density
# this fulfils the "regression" approach
# and allows exploration of how Ks(g) changes with temperature 

#NOTES:
#adding a term to fit N_init to all models. This feels really good, and predicted to observed looking better.
#Have done for 12_ and 13_ scripts.
#But model shapes are looking weird. 
#Also exploring allowing negative values of certain parameters, eg. d0 and Ks_slope, so_ks
#OK figure out weird-looking models next week!


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
IDE_kg<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks + N.Treatment)) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.18, b2=0.19, d0=-0.1, d1=0.15, d2=0.20, ks = 0.15), 
                lower=c(0,0,-5,0,0,-5,0), control = nls.control(maxiter=500, minFactor=1/204800000))))
#View(IDE_kg) # check that estimates aren't at boundaries.
write_csv(IDE_kg, "data-processed/IDE_kg")


#now take the root of this.

  newdata<-data.frame(Temperature=seq(3, 50, 0.1))
  
  for (j in 1:length(unique(IDE_kg$Species))){ #for each species
    this_species<-unique(IDE_kg$Species)[j]
    b1<-filter(IDE_kg, Species==this_species & term=="b1")$estimate
    b2<-filter(IDE_kg, Species==this_species & term=="b2")$estimate
    d0<-filter(IDE_kg, Species==this_species & term=="d0")$estimate
    d1<-filter(IDE_kg, Species==this_species & term=="d1")$estimate
    d2<-filter(IDE_kg, Species==this_species & term=="d2")$estimate
    ks<-filter(IDE_kg, Species==this_species & term=="ks")$estimate
    
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
  
  
  IDE_kg_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
  write_csv(IDE_kg_root, "data-processed/IDE_kg_root") #save
  
  #plot the root
  IDE_kg_root %>%
    ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
    facet_grid(~Species) +
    coord_cartesian(ylim = c(0, 10)) +
    theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(IDE_kg$Species))){ #for each species
  this_species<-unique(IDE_kg$Species)[j]
  b1<-filter(IDE_kg, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kg, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kg, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kg, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kg, Species==this_species & term=="d2")$estimate
  ks<-filter(IDE_kg, Species==this_species & term=="ks")$estimate
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=ks, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

IDE_umax_kg<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(IDE_umax_kg, "data-processed/IDE_umax_kg") #save

IDE_umax_kg<-read_csv("data-processed/IDE_umax_kg") #

#plot umax
IDE_umax_kg %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 2), xlim=c(0, 35)) +
  theme_bw() + ylab("Growth rate, r") 

#plot ks
IDE_umax_kg %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
IDE_kg_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                  ~ N_init + (((b1 * exp(b2*Temperature))
                               -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks + N.Treatment)) * day,
                  data= .,  start=list(N_init=5.7, b1 = 0.18, b2=0.19, d0=-0.1, d1=0.15, d2=0.20, ks = 0.15), 
                  lower=c(0,0,-5,0,0,-5,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                   ~ N_init + (((b1 * exp(b2*Temperature))
                                -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks + N.Treatment)) * day,
                   data= .,  start=list(N_init=5.7, b1 = 0.18, b2=0.19, d0=-0.1, d1=0.15, d2=0.20, ks = 0.15), 
                   lower=c(0,0,-5,0,0,-5,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 


#fitting IDE to cell density
#linear ks (=var)
IDE_kg_var<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.29, b2=0.12, d0=0.19, d1=0.13, d2=0.14, ks_int = 0, ks_slope=0.01), 
                lower=c(0,0,-5,0,0,-5,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))
View(IDE_kg_var)
write_csv(IDE_kg_var, "data-processed/IDE_kg_var")


#now take the root of this.
for (j in 1:length(unique(IDE_kg_var$Species))){ #for each species
  this_species<-unique(IDE_kg_var$Species)[j]
  b1<-filter(IDE_kg_var, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kg_var, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kg_var, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kg_var, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kg_var, Species==this_species & term=="d2")$estimate
  ks_int<-filter(IDE_kg_var, Species==this_species & term=="ks_int")$estimate
  ks_slope<-filter(IDE_kg_var, Species==this_species & term=="ks_slope")$estimate
  
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


IDE_kg_var_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(IDE_kg_var_root, "data-processed/IDE_kg_var_root") #save

#plot the root
IDE_kg_var_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(IDE_kg_var$Species))){ #for each species
  this_species<-unique(IDE_kg_var$Species)[j]
  b1<-filter(IDE_kg_var, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kg_var, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kg_var, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kg_var, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kg_var, Species==this_species & term=="d2")$estimate
  ks_int<-filter(IDE_kg_var, Species==this_species & term=="ks_int")$estimate
  ks_slope<-filter(IDE_kg_var, Species==this_species & term=="ks_slope")$estimate
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=ks_int+ks_slope*newdata$Temperature, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

umax_kg_var<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(umax_kg_var, "data-processed/umax_kg_var") #save

#plot umax
umax_kg_var %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
  theme_bw() + ylab("Predicted r-star, uM") 

#plot ks
umax_kg_var %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
IDE_kg_var_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                  ~ N_init + (((b1 * exp(b2*Temperature))
                               -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                  data= .,  start=list(N_init=5.7, b1 = 0.29, b2=0.12, d0=0.19, d1=0.13, d2=0.14, ks_int = 0, ks_slope=0.01), 
                  lower=c(0,0,-5,0,0,-5,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                   ~ N_init + (((b1 * exp(b2*Temperature))
                                -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + N.Treatment)) * day,
                   data= .,  start=list(N_init=5.7, b1 = 0.29, b2=0.12, d0=0.19, d1=0.13, d2=0.14, ks_int = 0, ks_slope=0.01), 
                   lower=c(0,0,-5,0,0,-5,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
  ggplot(aes(x=log.Particles.per.ml, y=.fitted, col=Species)) + geom_point() +
  facet_grid(~Species) +
  theme_bw() + ylab("observed") 

#
#
#
#maybe drop polynomial
#fitting IDE to cell density
#polynomial ks (=poly)
IDE_kg_poly<-alldata %>%
  group_by(Species) %>%
  do(tidy(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N.Treatment)) * day,
                data= .,  start=list(N_init=5.7, b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                lower=c(0,0,-5,0,0,-5,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))


write_csv(IDE_kg_poly, "data-processed/IDE_kg_poly")
View(IDE_kg_poly)

#now take the root of this.
for (j in 1:length(unique(IDE_kg_poly$Species))){ #for each species
  this_species<-unique(IDE_kg_poly$Species)[j]
  b1<-filter(IDE_kg_poly, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kg_poly, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kg_poly, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kg_poly, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kg_poly, Species==this_species & term=="d2")$estimate
  ks_int<-filter(IDE_kg_poly, Species==this_species & term=="ks_int")$estimate
  ks_slope<-filter(IDE_kg_poly, Species==this_species & term=="ks_slope")$estimate
  so_ks<-filter(IDE_kg_poly, Species==this_species & term=="so_ks")$estimate
  
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


IDE_kg_poly_root<-rbind(roots_AC, roots_CH, roots_CS, roots_TT) #bind together species data
write_csv(IDE_kg_poly_root, "data-processed/IDE_kg_poly_root") #save

#plot the root
IDE_kg_poly_root %>%
  ggplot(aes(x=temp, y=root, col=Species)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_bw() + ylab("Predicted r-star, uM") 

#now make plots of ks and umax across temp
for (j in 1:length(unique(IDE_kg_poly$Species))){ #for each species
  this_species<-unique(IDE_kg_poly$Species)[j]
  b1<-filter(IDE_kg_poly, Species==this_species & term=="b1")$estimate
  b2<-filter(IDE_kg_poly, Species==this_species & term=="b2")$estimate
  d0<-filter(IDE_kg_poly, Species==this_species & term=="d0")$estimate
  d1<-filter(IDE_kg_poly, Species==this_species & term=="d1")$estimate
  d2<-filter(IDE_kg_poly, Species==this_species & term=="d2")$estimate
  ks_int<-filter(IDE_kg_poly, Species==this_species & term=="ks_int")$estimate
  ks_slope<-filter(IDE_kg_poly, Species==this_species & term=="ks_slope")$estimate
  so_ks<-filter(IDE_kg_poly, Species==this_species & term=="so_ks")$estimate
  
  
  umax_ks<-data.frame(umax_temp=b1*exp(b2*newdata$Temperature)-(d0 + d1*exp(d2*newdata$Temperature)), 
                      ks_temp=ks_int+ks_slope*newdata$Temperature+so_ks*newdata$Temperature^2, Species=this_species, temp=newdata$Temperature)
  assign(paste("umax_ks",this_species,sep="_"), umax_ks)
}

IDE_umax_kg_poly<-rbind(umax_ks_AC, umax_ks_CH, umax_ks_CS, umax_ks_TT) #bind together species data
write_csv(IDE_umax_kg_poly, "data-processed/IDE_umax_kg_poly") #save

#plot umax
IDE_umax_kg_poly %>%
  ggplot(aes(x=temp, y=umax_temp, col=Species)) + geom_line() +
  coord_cartesian(ylim = c(0, 2)) +
  facet_grid(~Species) +
  theme_bw() + ylab("Predicted r-star, uM") 

#plot ks
IDE_umax_kg_poly %>%
  ggplot(aes(x=temp, y=ks_temp, col=Species)) + geom_line() +
  facet_grid(~Species) + theme_bw() + ylab("Predicted r-star, uM") 

#glance model summary
IDE_kg_poly_glance<-alldata %>%
  group_by(Species) %>%
  do(glance(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N.Treatment)) * day,
                data= .,  start=list(N_init=5.27, b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                lower=c(0,0,-5,0,0,-5,0,0,-1), control = nls.control(maxiter=500, minFactor=1/204800000))))

alldata %>%
  group_by(Species) %>%
  do(augment(nlsLM(log.Particles.per.ml 
                ~ N_init + (((b1 * exp(b2*Temperature))
                                   -(d0 + d1*exp(d2*Temperature)))) * (N.Treatment / (ks_int + ks_slope*Temperature + so_ks*Temperature^2 + N.Treatment)) * day,
                data= .,  start=list(N_init=5.27, b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, ks_int = 1, ks_slope=0.01, so_ks=0.5), 
                lower=c(0,0,-5,0,0,-5,0,0,-1), control = nls.control(maxiter=500, minFactor=1/204800000)))) %>%
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
