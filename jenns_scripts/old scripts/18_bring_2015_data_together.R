library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

TT_rstar<-cbind(read_csv("data-processed/TT_rstar.csv"), Species="TT")
CH_rstar<-cbind(read_csv("data-processed/CH_rstar.csv"), Species="CH")
AC_rstar<-cbind(read_csv("data-processed/AC_rstar.csv"), Species="AC")
CS_rstar<-cbind(read_csv("data-processed/CS_rstar.csv"), Species="CS")

all_rstar<-rbind(TT_rstar, CS_rstar, AC_rstar, CH_rstar)

TT_ks_umax<-cbind(read_csv("data-processed/TT_ks_umax.csv"), Species="TT")
CH_ks_umax<-cbind(read_csv("data-processed/CH_ks_umax.csv"), Species="CH")
AC_ks_umax<-cbind(read_csv("data-processed/AC_ks_umax.csv"), Species="AC")
CS_ks_umax<-cbind(read_csv("data-processed/CS_ks_umax.csv"), Species="CS")

all_ks_umax<-rbind(TT_ks_umax, CS_ks_umax, AC_ks_umax, CH_ks_umax)

3 10 17 24 31 38

all_ks_umax %>%
  #filter(term=="umax") %>%
  #group_by(term) %>%
  ggplot(aes(x = Temperature, y = estimate, color = factor(Species))) + geom_point(size = 4) +
  theme_bw() + geom_line() + facet_wrap(~term, scales = "free") + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)
ggsave("figures/all_ks_umax.png")

all_rstar %>%
  ggplot(aes(x = Temperature, y = median, color = factor(Species))) + geom_point(size = 4) +
  theme_bw() + geom_line() + labs(y="R-star") +
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.2)
ggsave("figures/all_rstar.png")
