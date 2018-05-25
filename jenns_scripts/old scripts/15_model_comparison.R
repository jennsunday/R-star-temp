#compare results from all models

#Compare AIC's of each model
all_model_glance<-rbind(
  mutate(IDE_kb_glance, model="IDE_kb"),
  mutate(IDE_kb_var_glance, model="IDE_kb_var"),
  mutate(IDE_kg_glance, model="IDE_kg"),
  mutate(IDE_kg_var_glance, model="IDE_kg_var"),
  mutate(NB_kg_glance, model="NB_kg"),
  mutate(NB_kg_var_glance, model="NB_kg_var"))

#View(all_model_glance)

#get the best models
all_model_glance %>%
  #filter(model %in% c("IDE_kb", "IDE_kb_var", "IDE_kg", "IDE_kg_var")) %>%
  group_by(Species)  %>%
  summarize(min_AIC=min(AIC), min_AIC_model=model[which(AIC == min(AIC))])
  #summarize(max_lik=max(logLik), max_lik_model=model[which(logLik == max(logLik))])

write_csv(all_model_glance, "data-processed/all_model_glance.csv")

#Read in all roots
all_model_roots<-rbind(
  mutate(read_csv("data-processed/IDE_kg_root"), model="IDE", resource_affects="growth", ks="fixed"),
  mutate(read_csv("data-processed/IDE_kg_var_root"), model="IDE", resource_affects="growth", ks="var"),
  mutate(read_csv("data-processed/IDE_kb_root"), model="IDE", resource_affects="birth", ks="fixed"),
  mutate(read_csv("data-processed/IDE_kb_var_root"), model="IDE", resource_affects="birth", ks="var"),
  mutate(read_csv("data-processed/NB_kg_root"), model="NB", resource_affects="growth", ks="fixed"),
  mutate(read_csv("data-processed/NB_kg_var_root"), model="NB", resource_affects="growth", ks="var"))

#plot all the roots
all_model_roots %>%
  #filter(ks=="var") %>%
  ggplot(aes(x=temp, y=root, colour=model, linetype=ks)) + geom_line() +
  coord_cartesian(ylim = c(0, 10)) +
  facet_grid(resource_affects~Species) +
  theme_bw() + ylab("Predicted r-star, uM") 
ggsave("figures/all_models_roots.pdf", width=8, height=4)


#plot the best roots
  all_model_roots %>%
  filter(model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="AC" |
          model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="CH"|
           model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="CS"|
            model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="TT"|
           model=="NB" & resource_affects=="growth" & ks=="var" & Species=="TT") %>%
  ggplot(aes(x=temp, y=root, colour=Species, linetype=ks)) + geom_line() + 
    facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 10), xlim=c(-1,50))  +
  theme_bw() + ylab("R-star, uM") 
ggsave("figures/all_models_roots.pdf", width=8, height=2)


#plot umax
umax_NB_kg_var<-read_csv("data-processed/umax_NB_kg_var") #
umax_NB_kg<-read_csv("data-processed/umax_NB_kg") #
umax_NB_all<-rbind(mutate(umax_NB_kg_var, ks="var"), 
                   mutate(umax_NB_kg, ks="fixed"))

umax_NB_all %>%
  filter(ks=="fixed" & Species=="AC" |
         ks=="fixed" & Species=="CH"|
         ks=="fixed" & Species=="CS"|
         ks=="fixed" & Species=="TT"|
         ks=="var" & Species=="TT") %>%
  ggplot(aes(x=temp, y=umax_temp, colour=Species, linetype=ks)) + geom_line() + 
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 1.6), xlim=c(-1,50))  +
  theme_bw() + ylab("umax") 
ggsave("figures/all_models_umax.pdf", width=8, height=2)

#plot predicted and observed nutrient draw-downn
Nitrate_orig<-4.362499 # N original from R* experiment 2

# Read in N equilibrium data from R* experiment 2
rstar_by_temp<-read_csv("data-processed/rstar_by_temp.csv") %>% 
  filter(term=="a") %>%   
  #filter(!temp==38 | !Species=="TT")%>%
  #filter(!temp==38 | !Species=="CS")%>%
  #filter(!temp==38 | !Species=="AC")%>%
  mutate(draw_down=Nitrate_orig-estimate)

all_model_roots<-all_model_roots %>%
  mutate(N_drawdown=ifelse(root>Nitrate_orig, 0, Nitrate_orig-root)) %>%
  filter(model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="AC" |
         model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="CH"|
        model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="CS"|
        model=="NB" & resource_affects=="growth" & ks=="fixed" & Species=="TT"|
        model=="NB" & resource_affects=="growth" & ks=="var" & Species=="TT")

all_model_roots%>%
  ggplot(aes(x=temp, y=root, colour=model, linetype=ks, model=resource_affects)) + geom_line() +
  coord_cartesian(ylim = c(0, 5), xlim=c(3,50)) +
  facet_grid(~Species) +
  theme_bw() + ylab("Available N at equilibrium, uM") + 
  geom_point(data=rstar_by_temp, aes(x = temp, y = estimate), inherit.aes=F)


rstar_by_temp %>% 
  ggplot(aes(x = temp, y = draw_down, color=as.factor(Species))) + geom_point(size = 2) +
  theme_bw() + facet_grid(~Species) + geom_smooth(method = 'loess') + ylab("Nitrate draw down at equilibrium, uM")
ggsave("figures/rstar_with_temp_2017.pdf", width=7, height=2)


