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


all_fits<-rbind(read_csv("data-processed/all_summary_fit_IDE_kb_fixed.csv") %>% mutate(model="IDE_kb_fixed"),
                read_csv("data-processed/all_summary_fit_IDE_kg_fixed.csv") %>% mutate(model="IDE_kg_fixed"),
                read_csv("data-processed/all_summary_fit_IDE_kb_var.csv") %>% mutate(model="IDE_kb_var"),
                read_csv("data-processed/all_summary_fit_IDE_kg_var.csv") %>% mutate(model="IDE_kg_var"),
                read_csv("data-processed/all_summary_fit_NB_fixed.csv") %>% mutate(model="NB_fixed"),
                read_csv("data-processed/all_summary_fit_NB_var.csv") %>% mutate(model="NB_var"))
chart<-all_fits %>%
  select(Species, model, AIC, logLik)
write.csv(chart, "data-processed/model_AICs.csv")

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

#write out results
write_csv(best_model_params, "data-processed/best_model_params.csv")
write_csv(best_model_params_booted, "data-processed/best_model_params_booted.csv")

p1<-ggplot() + geom_line(data=best_model_params, aes(x=Temperature, y=root, colour=Species)) + 
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 7), xlim=c(-1,50))  +
  theme_bw() + ylab("R-star, uM") +
  geom_ribbon(data=best_model_params_booted, 
              aes(x=Temperature, ymin = root_lwr_CI, ymax = root_upr_CI), alpha = .1)

#read in partial-model results
indirect_umax<-read_csv("data-processed/indirect_umax_ks.csv") %>% 
  filter(term=="umax")

p2<-ggplot() + geom_line(data=best_model_params, aes(x=Temperature, y=umax, colour=Species)) + 
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 1.9), xlim=c(-1,50))  +
  theme_bw() + ylab("umax, day-1") +
  geom_ribbon(data=best_model_params_booted, 
              aes(x=Temperature, ymin = umax_lwr_CI, ymax = umax_upr_CI), alpha = .1) +
  geom_point(data=indirect_umax, aes(x=Temperature, y=estimate, colour=Species)) +
  geom_errorbar(data=indirect_umax, aes(x=Temperature, ymin=estimate-std.error, 
                                        ymax=estimate+std.error, colour=Species, width=0))

together<-grid.arrange(p1, p2, nrow = 2)
ggsave("figures/best_roots_and_umax.pdf", together, width=8, height=4)


#goal: add some info to this plot that shows the rate of minimum R* compared to the range of maximum growth rate
#first: where is the minimum R* / maximum r?
head(best_model_params)

minmax<-best_model_params %>%
    group_by(Species)  %>%
    summarize(min_root=min(root), Tminroot=Temperature[which(root==min(root))],
              max_umax=max(umax), Tmaxumax=Temperature[which(umax==max(umax))])
lowroots$Species

lowroots<-best_model_params %>%
  group_by(Species)  %>%
  mutate(lowroot=min(root)+((5-min(root))/100)) %>%
  filter(root<lowroot) %>%
  mutate(ydum=1)

highrs<-best_model_params %>%
  filter(umax>0) %>%
  group_by(Species)  %>%
  mutate(highr=max(umax)-(max(umax)/100)) %>%
  filter(umax>highr) %>%
  mutate(ydum=0)

p3<-ggplot() + geom_line(data=highrs, aes(y=ydum, x=Temperature, col=Species, size=0.8)) +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(-1, 2), xlim=c(-1,50))  +
  geom_line(data=lowroots, aes(y=ydum, x=Temperature, col=Species, size=0.8)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
  

together<-grid.arrange(p1,p3, p2, nrow = 3)
ggsave("figures/best_roots_and_umax.pdf", together, width=8, height=7)


#Fig. 3 display fits of umax of  best-fit models to indirect parameter estimates from other scripts#####
predict_growth_NB<-rbind(read_csv("data-processed/predict_growth_NB_fixed.csv") %>% 
                           select(N, Temp, Species, r_pred) %>%
                           mutate(ks_type="fixed"),
                         read_csv("data-processed/predict_growth_NB_var.csv") %>% 
                           select(N, Temp, Species, r_pred) %>%
                           mutate(ks_type="var"))
indirect_r_fit<-rbind(read_csv("data-processed/indirect_r.csv") %>% 
                        mutate(ks_type="fixed"),
                         read_csv("data-processed/indirect_r_NB_var_compare.csv") %>% 
                        mutate(ks_type="var"))


growth_predicted_with_best_models<-predict_growth_NB %>%
  filter(ks_type=="var" & Species=="AC" |
           ks_type=="fixed" & Species=="CH"|
           ks_type=="fixed" & Species=="CS"|
           ks_type=="var" & Species=="TT")

indirect_r_fit_best<-indirect_r_fit %>%
  filter(ks_type=="var" & Species=="AC" |
           ks_type=="fixed" & Species=="CH"|
           ks_type=="fixed" & Species=="CS"|
           ks_type=="var" & Species=="TT")


Topt_best_model<-growth_predicted_with_best_models %>%
  group_by(Species,N) %>%
  summarize(max_r=max(r_pred), Topt=Temp[which(r_pred==max(r_pred))])

P1<-indirect_r_fit_best %>%
  ggplot(aes(y=estimate, x=as.numeric(N.Treatment), col=as.factor(Temperature))) +
  facet_grid(~Species, scales="free_y") + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))+
  geom_line(data=filter(growth_predicted_with_best_models, Temp %in% unique(alldata$Temperature)), aes(y=r_pred, x=N, colour=as.factor(Temp))) +
  coord_cartesian(ylim = c(0, 1.9))  +
  theme_bw()

P2<-indirect_r_fit_best %>%
  ggplot(aes(y=estimate, x=Temperature, col=as.factor(N.Treatment))) +
  facet_grid(~Species, scales="free_y") + geom_point() +
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error, width=0))+
  geom_line(data=growth_predicted_with_best_models, aes(y=r_pred, x=Temp, colour=as.factor(N))) +
  geom_point(data=Topt_best_model, aes(y=max_r, x=Topt, colour=as.factor(N)), inherit.aes=FALSE, shape=21, fill=grey(0.3)) +
  coord_cartesian(ylim = c(0, 1.9), xlim=c(13,28))  +
  theme_bw()

together<-grid.arrange(P1, P2, nrow = 2)
ggsave("figures/best_growth_fits.pdf", together, width=8, height=4)
