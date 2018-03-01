#libraries
library(bbmle)
library(tidyverse)
library(stringr)

#hypothetical_curve_1<-cbind(read_csv("data-processed/Norberg_fits_by_N_TT_predictions.csv"), Species="TT") %>%
#  filter(N.Treatment==440)

hypothetical_curve_1<-cbind(read_csv("data-processed/Norberg_fits_by_N_CH_predictions.csv"), Species="1") %>%
  filter(N.Treatment==440) %>%
  mutate(r_predict=r_predict*0.6)

hypothetical_curve_2<-cbind(read_csv("data-processed/Norberg_fits_by_N_CH_predictions.csv"), Species="2") %>%
  filter(N.Treatment==440) %>%
  mutate(temperature=temperature+6)

head(together)
together<-rbind(hypothetical_curve_1, hypothetical_curve_2)

k.s<-3
m<-0.1
together<-together %>%
  mutate(r_star=m*k.s/(r_predict-m))#R* = mKs / umax - m

together %>%
  filter(r_predict>0) %>%
  ggplot(aes(y=r_predict, x=temperature, color=Species)) + 
  coord_cartesian(ylim = c(0, 8)) + geom_line(size=1) + 
  theme_bw()
ggsave("figures/hypothetical_TPC.pdf", width=3, height=2)

together %>%
  filter(r_predict>0.1) %>%
  ggplot(aes(y=r_star, x=temperature, color=Species)) + 
  coord_cartesian(ylim = c(0, 8)) + geom_line(size=1) + 
  theme_bw()
ggsave("figures/hypothetical_r_star.pdf", width=3, height=2)

