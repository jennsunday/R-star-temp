#goal: 
# fit a monod curve and growth rate simultaneously to all of the cell density data in each temperature
# this fulfils the "regression" approach
# and allows exploration of how Ks changes with temperature 
# (but does not simultaneously fit growth rates across temperatures because don't know how Ks changes with temp)

library(broom)
library(purrr)
library(tidyverse)
library(minpack.lm)
library(nlstools)

TTfilteredN<-read_csv("data-processed/TTfilteredN.csv")

#fitting monod to cell density
test<-TTfilteredN %>%
    group_by(Temperature) %>%
    do(tidy(nlsLM(log.Particles.per.ml ~ 6 + umax * (N.Treatment / (ks+ N.Treatment))*day,
                  data= .,  start=list(ks = 6.6, umax = 0.59), 
                  lower=c(0.001, 0))))

tempdata<-TTfilteredN %>%
  filter(Temperature==unique(TTfilteredN$Temperature)[1])
testboot<-nlsBoot(nlsLM(log.Particles.per.ml ~ 6 + umax * (N.Treatment / (ks+ N.Treatment))*day,
                        data= tempdata,  start=list(ks = 6.6, umax = 0.59), algorithm="port", 
                        lower=c(0.001, 0)), niter=999)
summary(testboot)

test<-nlsLM(log.Particles.per.ml ~ 6 + umax * (N.Treatment / (ks+ N.Treatment))*day,
        data= TTfilteredN,  start=list(ks = 6.6, umax = 0.59), algorithm="port", 
        lower=c(0.001, 0))

nlsBoot(test,  niter=999)

test_pred<-TTfilteredN %>%
  group_by(Temperature) %>%
  do(augment(nlsLM(log.Particles.per.ml ~ 6 + umax * (N.Treatment / (ks+ N.Treatment))*day,
                data= .,  start=list(ks = 6.6, umax = 0.59), algorithm="port", 
                lower=c(0.001, 0),
                control = nls.control(maxiter=500, minFactor=1/204800000)), conf.int=T))

test_pred %>%
  ggplot(aes(y=log.Particles.per.ml, x=day)) + facet_grid(Temperature~N.Treatment) +
  geom_point() + geom_line(aes(y=.fitted, color="red"))  + theme_bw() +
  #geom_ribbon(aes(ymax=.fitted+.resid, ymin=.fitted-.resid))
ggsave("figures/fit_monod_to_growth_data.pdf")

test %>%
  ggplot(aes(y=estimate, x=Temperature)) + geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) +
  facet_wrap(~term, scales="free")
ggsave("figures/TT_ks_umax_fit_to_cell_dens.pdf", widge)

?nlsBoot

#define monod curve function for growth rate
monodcurve<-function(N.Treatment, ks, umax){
  growth_rate<-umax * (N.Treatment / (ks+ N.Treatment))
  growth_rate}

#define growth ratefunction for log cell density
growthcurve<-function(day, umax, ks, int, N.Treatment){
  log.cell.dens<- int + umax * (N.Treatment / (ks+ N.Treatment))*day
  log.cell.dens}

#create dummy data for displaying monod curve in ggplot
newdataN=seq(0, 440, 1)



test %>%
  filter(term=="ks") %>%
  ggplot(aes(y=estimate, x=Temperature)) + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)


test2<-data.frame()
 pred=monodcurve(newdataN, test$estimate[1], test$estimate[2])) %>%
  gglot(aes=)


#create dummy data for displaying cell density curve in ggplot
pred_growth_df<-expand.grid(day=seq(0, 6, 1), N.Treatment=c(0, 11, 22, 55, 110, 220, 330, 440),
                            temperature=c(13, 16, 19, 22, 25, 28)) 

pred_growth_df %>%
  mutate(pred_growth=growthcurve(newdata$day, 
                          test$estimate[1],
                          test$estimate[2], test$estimate[3], newdata$N.Treatment))


