library(broom)
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(tidyverse)
library(nls.multstart)
library(gridExtra)

#compare model fits
all_fits<-rbind(read_csv("data-processed/all_summary_fit_IDE_kb_fixed.csv") %>% mutate(model="IDE_kb_fixed"),
                read_csv("data-processed/all_summary_fit_IDE_kg_fixed.csv") %>% mutate(model="IDE_kg_fixed"),
                read_csv("data-processed/all_summary_fit_IDE_kb_var.csv") %>% mutate(model="IDE_kb_var"),
                read_csv("data-processed/all_summary_fit_IDE_kg_var.csv") %>% mutate(model="IDE_kg_var"),
                read_csv("data-processed/all_summary_fit_NB_fixed.csv") %>% mutate(model="NB_fixed"),
                read_csv("data-processed/all_summary_fit_NB_var.csv") %>% mutate(model="NB_var"))

all_fits %>%
  group_by(Species)  %>%
  summarize(min_AIC=min(AIC), min_AIC_model=model[which(AIC == min(AIC))])

#Read in NB parameters
all_model_roots<-rbind(
  mutate(read_csv("data-processed/all_roots_ks_umax_NB_fixed.csv"), model="NB", resource_affects="growth", ks="fixed"),
  mutate(read_csv("data-processed/all_roots_ks_umax_NB_var.csv"), model="NB", resource_affects="growth", ks="var"))

#Read in NB booted parameters
all_model_boot_confint<-rbind(
  mutate(read_csv("data-processed/all_boots_confint_NB_fixed"), model="NB", resource_affects="growth", ks="fixed"),
  mutate(read_csv("data-processed/all_boots_confint_NB_var"), model="NB", resource_affects="growth", ks="var"))

#just the best-fit models
best_model_params<-all_model_roots %>%
  filter(ks=="var" & Species=="AC" |
           ks=="fixed" & Species=="CH"|
           ks=="fixed" & Species=="CS"|
           ks=="var" & Species=="TT")

best_model_params_booted<-all_model_boot_confint %>%
  filter(ks=="var" & Species=="AC" |
           ks=="fixed" & Species=="CH"|
           ks=="fixed" & Species=="CS"|
           ks=="var" & Species=="TT")

#read in partial-model results
indirect_umax<-read_csv("data-processed/indirect_umax_ks.csv") %>% 
  filter(term=="umax")

p1<-ggplot() + geom_line(data=best_model_params, aes(x=Temperature, y=root, colour=Species)) + 
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 10), xlim=c(-1,50))  +
  theme_bw() + ylab("R-star, uM") +
  geom_ribbon(data=best_model_params_booted, 
              aes(x=Temperature, ymin = root_lwr_CI, ymax = root_upr_CI), alpha = .1)

p2<-ggplot() + geom_line(data=best_model_params, aes(x=Temperature, y=umax, colour=Species)) + 
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 1.5), xlim=c(-1,50))  +
  theme_bw() + ylab("umax, day-1") +
  geom_ribbon(data=best_model_params_booted, 
              aes(x=Temperature, ymin = umax_lwr_CI, ymax = umax_upr_CI), alpha = .1) +
  geom_point(data=indirect_umax, aes(x=Temperature, y=estimate, colour=Species)) +
  geom_errorbar(data=indirect_umax, aes(x=Temperature, ymin=estimate-std.error, 
                                        ymax=estimate+std.error, colour=Species, width=0))

together<-grid.arrange(p1, p2, nrow = 2)
ggsave("figures/best_roots_and_umax.pdf", together, width=8, height=4)
