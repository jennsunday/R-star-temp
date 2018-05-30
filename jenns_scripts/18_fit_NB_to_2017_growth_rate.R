cell_results<-read.csv("data-processed/logistic_growth_fits_rep_r-star.csv") %>%
  mutate(r_corr=r+0.1) #correct for an added death rate of 0.1 per day

cell_results_within<-cell_results %>%
  filter(Species %in% c("AC", "CS", "TT") & temp!=38 | 
           Species %in% c("CH"))
#fitting NB to growth rate
#fixed k
####################################
head(cell_results)

all_params_fit<-data.frame()
all_summary_fit<-data.frame()
all_preds_fit<-data.frame()
all_roots_ks_umax<-data.frame()

for(i in unique(cell_results$Species)){
  data<-filter(cell_results_within, Species==i)
  fit<-nls_multstart(r_corr ~ a*exp(b*temp)*(1-((temp-z)/(w/2))^2),
                     data= data,  iter = 250,
                     start_lower = c(a = 0.1, b=0.0001, z=5, w=5),
                     start_upper = c(a = 0.6, b=0.2, z=40, w=40),
                     supp_errors = 'Y',
                     convergence_count = 100,
                     na.action = na.omit,
                     lower=c(a = 0, b=0, z=-50, w=0))
  
  params_fit<-tidy(fit) %>% mutate(Species=i)
  all_params_fit<-rbind(all_params_fit, params_fit)
  
  summary_fit<-glance(fit) %>% mutate(Species=i)
  all_summary_fit<-rbind(all_summary_fit, summary_fit)
  
  preds_fit<- augment(fit)  %>% mutate(Species=i)
  all_preds_fit<-rbind(all_preds_fit, preds_fit)
}

#predict results

predict_grid<-expand.grid(Species=unique(cell_results$Species), Temperature=newdata$Temperature)

backtogether<-data.frame()
for(i in 1:4){
  predict_grid_sp<-predict_grid %>%
    filter(Species==unique(alldata$Species)[i]) 
  NB_fit_sp<-all_params_fit %>%
    filter(Species==unique(alldata$Species)[i]) 
  predict_grid_sp<-predict_grid_sp%>% 
    mutate(a=filter(NB_fit_sp, term=="a")$estimate,
           b=filter(NB_fit_sp, term=="b")$estimate,
           z=filter(NB_fit_sp, term=="z")$estimate,
           w=filter(NB_fit_sp, term=="w")$estimate)
  backtogether<-rbind(backtogether, predict_grid_sp)
}
predict_grid_full<-backtogether


#write function
growth_function_17<-function(Temp,a,b,z,w){
  a*exp(b*Temp)*(1-((Temp-z)/(w/2))^2)}


predict_growth_17 <- predict_grid_full %>%
  mutate(r_pred=growth_function_17(predict_grid_full$Temp,
                                predict_grid_full$a,
                                predict_grid_full$b,
                                predict_grid_full$z,
                                predict_grid_full$w))

#bring in direct estimates of R-star for multi-panel display
write_csv(rstar_by_temp, "data-processed/rstar_by_temp.csv")

rstar_by_temp<-read_csv("data-processed/rstar_by_temp.csv") %>% 
  filter(term=="a")


P2<-ggplot() + geom_line(data=predict_growth_17, aes(y=r_pred, x=Temperature, col=Species)) +
  facet_grid(~Species) +
  geom_point(data=cell_results, aes(y=r_corr, x=temp, col=Species)) +
  coord_cartesian(ylim = c(-.2, 2.2), xlim = c(0, 50)) + 
  ylab("growth rate, d-1") +
  theme_bw()


P1<-rstar_by_temp %>% 
  filter(!temp==38 | !Species=="TT")%>%
  filter(!temp==38 | !Species=="CS")%>%
  filter(!temp==38 | !Species=="AC")%>%
  ggplot(aes(x = temp, y = estimate, color=Species)) + geom_point(size = 1) +
  theme_bw() + facet_grid(~Species) + geom_smooth(method = 'loess') + ylab("R-star, uM") +
  coord_cartesian(xlim = c(0, 50)) + 
  xlab("Temperature")

together<-grid.arrange(P1, P2, nrow = 2)
ggsave("figures/2017_direct_estimates.pdf", together, width=8, height=4)

