#plan: 
# fit a monod curve to all of the cell density data in each temperature - this fulfils the "regression" approach
# and allows exploration of how Ks changes with temperature (but does not connect growth rates across temperatures)
# fit a double-exponential TPC model to umax to extrapolate
# fit a GAM to ks to extrapolate (reread Bestian) (but be upfront that 2015 data were not great for estimating Ks)
# then consider using the big model, all parameters b1, b2, d0, d1, d2, ks=f(temp), in order to 
# 1) plot the 2-d figure from Thomas for each species.
# 2) use uniroot to plot R* ~ temp for each species.


# for Nitrate, refit decay model using Michaelis Menton
# think hard/read about whether equilibrial environmental nitrate in continuous flow should be higher than 
# R*. E.g. if R* can be calcated as the root of Thomas equation, where is the mortality component of the continuous flow?
# hmmm, ok, I thought:
# I think the mortality component is important because R*


library(broom)
library(purrr)
library(tidyverse)
library(minpack.lm)
library(nlstools)

#Joey's code fitting Norberg to cell density
nlsLM(cell_density ~ 800 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
      data= cd,  
      start=list(z=df$z[[1]],w=35,a= df$a[[1]], b=df$b[[1]]),
      lower = c(0, 0, -0.2, -0.2),
      upper = c(30, 80, 5, 0.5),
      control = nls.control(maxiter=1024, minFactor=1/204800000))


#Jenn's old code fitting monod to growth rate
nlsLM(r_estimate ~ umax* (N.Treatment / (ks+ N.Treatment)),
            data= .,  start=list(ks = 1, umax = 1), algorithm="port", lower=list(c=0.01, d=0),
            control = nls.control(maxiter=500, minFactor=1/204800000)

#fitting monod to cell density - just TT, just 13
for(i in 1:100){
test<-TTfilteredN %>%
  filter(Temperature==13) %>%
  #group_by(N.Treatment) %>% 
  #sample_n(., 6, replace=FALSE) %>%
  do(tidy(nlsLM(log.Particles.per.ml ~ int + umax * (N.Treatment / (ks+ N.Treatment))*day,
    data= .,  start=list(ks = 6.6, umax = 0.59, int=6.7), algorithm="port", 
    lower=c(0.001, 0, 2),
    control = nls.control(maxiter=500, minFactor=1/204800000))))
  temp<-data.frame(slope=test[2,2],
                 intercept=onefit[1,2],
                 N.Treatment=unique(TTfilteredN$N.Treatment)[j],
                 Temperature=unique(TTfilteredN$Temperature)[k])
}

#reminder of the function for cell number
log_cell_number = int+growthrateslope*day
#reminder of the fuction for growth rate
growthrate = (b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))

umaxcurve<-function(Temperature, b1,b2,d0, d1, d2){
  growth_rate<-(b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))
  growth_rate}

celldenscurve<-function(Temperature, int, b1,b2,d0, d1, d2){
  growth_rate<- int + ((b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))) * 7
  growth_rate}

TTfilteredN %>%
  filter(N.Treatment==440) %>%
  plot(log.Particles.per.ml~Temperature, data=., ylim=c(0,12))
curve(celldenscurve(x,5.5, 0.29, 0.112, 0, 0.135, 0.135),col='red', lwd=2, xlim=c(1, 38), add=T)

#int = 5.5
#b1 0.29
#b2 0.112
#d0 0
#d1 0.135
#d2 0.135

celldenscurve<-function(Temperature, int, b1,b2,d0, d1, d2,day){
  growth_rate<- int + ((b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))) * day
  growth_rate}

for(i in 1:100){
  test_de_220<-TTfilteredN %>%
    filter(N.Treatment==220) %>%
    do(tidy(nlsLM(log.Particles.per.ml ~ int + 
                    ((b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))) * day,
                  data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135, int=5.5), algorithm="port", 
                  lower=c(0,0,0,0,0,3), control = nls.control(maxiter=500, minFactor=1/204800000))))
  temp<-data.frame(slope=test[2,2],
                   intercept=onefit[1,2],
                   N.Treatment=unique(TTfilteredN$N.Treatment)[j],
                   Temperature=unique(TTfilteredN$Temperature)[k])
}



TTfilteredN %>%
  filter(N.Treatment==440) %>%
  plot(log.Particles.per.ml~Temperature, data=., ylim=c(0,20))

thiscurve<-test_de_440
curve(celldenscurve(x, b1=thiscurve$estimate[1], 
                      b2=thiscurve$estimate[2], 
                      d0=thiscurve$estimate[3], 
                      d1=thiscurve$estimate[4], 
                      d2=thiscurve$estimate[5], 
                      int=thiscurve$estimate[6], day=5),
      col='red', lwd=2, xlim=c(1, 38), add=T)

curve(umaxcurve(x, b1=test_de[1,2], b2=test_de[2,2], d0=test_de[3,2], d1=test_de[4,2], d2=test_de[5,2]),
      col='red', lwd=2, xlim=c(1, 38), ylim=c(-0.5,2))

head(TTfilteredN)

    N(t) = N0e^rt
    
nbcurve<-function(temp,z,w,a,b){
  res<-a*exp(b*temp)*(1-((temp-z)/(w/2))^2) * 
  res
}

fit_growth <- function(df){
  res <- try(nlsLM(cell_density ~ 800 * exp((a*exp(b*temp)*(1-((temp-z)/(w/2))^2))*(days)),
                   data= cd,  
                   start=list(z=df$z[[1]],w=df$w[[1]],a= df$a[[1]], b=df$b[[1]]),
                   lower = c(0, 0, -0.2, -0.2),
                   upper = c(30, 80, 1.2, 0.3),
                   control = nls.control(maxiter=1024, minFactor=1/204800000)))
  if(class(res)!="try-error"){
    out1 <- tidy(res) %>% 
      select(estimate, term) %>% 
      spread(key = term, value = estimate)
    out2 <- glance(res)
  }
  all <- bind_cols(out1, out2)
  all
}

b1vals<-seq(-3,3,0.5)
b2vals<-seq(0.02,1.2,0.3)
Ksvals<-seq(0.01, 8, 1)
d0vals<-seq(0.02,1.2,0.3)
d1vals<-seq(0.01,5,0.5)
d2vals<-seq(0.02,1.2,0.2)

df <-  expand.grid(b1 = b1vals, b2 = b2vals, Ks = Ksvals, d0 = d0vals, d1 =d1vals, d2=d2vals) %>% 
  mutate(unique_id = rownames(.))
  
fit_growth <- function(df){
temp_nitrate_model<-try(nlsLM(growth_rate 
      ~ (b1 * exp(b2*Temperature) * N.Treatment/(N.Treatment+Ks))-(d0 + d1*exp(d2*Temperature)),
      data=TT_data, 
      start=list(b1=df$b1[[1]], b2=df$b2[[1]], Ks=df$Ks[[1]], d0=df$d0[[1]], d1=df$d1[[1]], d2=df$d2[[1]]),
      upper = c(5, 5, 12, 0, 0, 0),
      control = nls.control(maxiter=1024, minFactor=1/204800000)))
if(class(temp_nitrate_model)!="try-error"){
  out1 <- tidy(temp_nitrate_model) %>% 
    select(estimate, term) %>% 
    spread(key = term, value = estimate)
  out2 <- glance(temp_nitrate_model)
}}

all <- bind_cols(out1, out2)
all




df_split <- df %>% 
  split(.$unique_id)

output <- df_split %>%
  map_df(fit_growth, .id = "run") 

return(output)



samples <- rep(1, 10) 

bootstrap_time_series_fits <- samples %>% ### this step gets us the replication
  map_df(resample, .id = "replicate")


coef(temp_nitrate_model)



umaxcurve<-function(Temperature, b1,b2,d0, d1, d2){
  growth_rate<-(b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))
  growth_rate}

-0.1*exp(30*-0.05)
curve(umaxcurve(x,0.12, 0.06, 0, 0.055, 0.08),col='red', lwd=2, xlim=c(1, 38))
abline(0,0)
#b1 0.12, 0.1, 0.09
#d1 0.035, 0.04, 0.05, 0.06
#d0 - lowers height of curve
#b2 - 0.065, 0.06


  
(b1 * exp(b2*Temperature) * N.Treatment/(N.Treatment+Ks))-(d0 + d1*exp(d2*Temperature))

umaxcurve(coef(temp_nitrate_model)[1], coef(temp_nitrate_model)[2], coef(temp_nitrate_model)[3], 
          coef(temp_nitrate_model)[4], coef(temp_nitrate_model)[5])

minpack.lm

nls.tools

TT_data<-lmfits %>%
  filter(term=="day") %>%
  rename(growth_rate=estimate)

names(TT_data)
