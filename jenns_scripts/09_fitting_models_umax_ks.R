#goal: fit DDE to umax, fit GAM to ks
library(broom)
library(purrr)
library(tidyverse)
library(minpack.lm)
library(nlstools)
#don't let ks go below zero
library(mgcv)

monod_growth_fit<-read_csv("data-processed/monod_growth_fit.csv")

#DDE model to umax
#get parameters
test_dde<-monod_growth_fit %>%
    filter(term=="umax") %>%
    group_by(Species) %>%
    do(tidy(nlsLM(estimate ~ 
                    ((b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))),
                  data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135), algorithm="port", 
                  lower=c(0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000))))

#plot lines
newdata<-data.frame(Temperature=seq(3, 38, 0.5))
test_dde_pred<-monod_growth_fit %>%
  filter(term=="umax") %>%
  group_by(Species) %>%
  do(augment(nlsLM(estimate ~ 
                  ((b1 * exp(b2*Temperature))-(d0 + d1*exp(d2*Temperature))),
                data= .,  start=list(b1 = 0.29, b2=0.112, d0=0, d1=0.135, d2=0.135), algorithm="port", 
                lower=c(0,0,0,0,0), control = nls.control(maxiter=500, minFactor=1/204800000)), newdata=newdata))

monod_growth_fit %>%
  filter(term=="umax") %>%
  ggplot(aes(y=estimate, x=Temperature, col=Species)) + geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) +
  facet_grid(~Species, scales="free_y") +
  geom_line(data=test_dde_pred, aes(y=.fitted)) +
  coord_cartesian(ylim = c(0, 2))
ggsave("figures/umax_fit.pdf", width=7, height=2)


#GAM model to ks
test_gam<-monod_growth_fit %>%
  filter(term=="ks") %>%
  group_by(Species) %>%
  do(tidy(gam(estimate~s(Temperature, bs="cr", k=5), data=.)))

test_gam_pred<-monod_growth_fit %>%
  filter(term=="ks") %>%
  group_by(Species) %>%
  do(augment(gam(estimate~s(Temperature, bs="cr", k=5), data=.), newdata=newdata))

monod_growth_fit %>%
  filter(term=="ks") %>%
  ggplot(aes(y=estimate, x=Temperature, col=Species)) + geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2) +
  facet_grid(~Species, scales="free_y") +
  geom_line(data=test_gam_pred, aes(y=.fitted)) +
  geom_ribbon(data=test_gam_pred, aes(y=.fitted, ymax=(.fitted+.se.fit), ymin=(.fitted-.se.fit)),
              alpha=0.2, linetype=0)
ggsave("figures/ks_fit.pdf", width=7, height=2)

