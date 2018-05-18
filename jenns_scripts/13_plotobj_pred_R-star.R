
rstar_by_temp<-read_csv("data-processed/rstar_by_temp.csv")
allspecies_roots<-read_csv("data-processed/allspecies_roots.csv")

rstar_by_temp_sub<-rstar_by_temp %>% 
  filter(term=="a") %>% 
  filter(!temp==38 | !Species=="TT")%>%
  filter(!temp==38 | !Species=="CS")%>%
  filter(!temp==38 | !Species=="AC")

  ggplot(aes(x = temp, y = estimate, color=as.factor(Species))) + geom_point(size = 2) +
  theme_bw() + facet_grid(~Species) + geom_smooth(method = 'loess') + ylab("Nitrate at equilibrium, uM")

allspecies_roots %>%
  ggplot(aes(x=temp, y=root, col=model)) + geom_line() +
  facet_grid(~Species) +
  coord_cartesian(ylim = c(0, 10)) + 
  theme_bw() + ylab("Predicted r-star, uM") +
  geom_point(data=rstar_by_temp_sub, aes(x = temp, y = estimate), inherit.aes=F, size = 2)
ggsave("figures/Predicted r-star.pdf", width=7, height=2)