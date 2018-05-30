#Goal: create functions to fit logistic growth curve to nitrate data

#read in libraries
library(simecol)
library(tidyverse)
library(broom)
library(minpack.lm)
library(grid)

#read in data
cdata<-read_csv("data-processed/processed_nitrate_data_SE.csv")
names(cdata)

#remove a rep with missing data
cdata<- cdata %>% 
  filter(temp!=10 | rep!=4 | species!=2)

#figure out the mean starting condition
N_init<-cdata$n.fixed[1]

#define N
#note: "n.fixed" is nitrate column when day 0 N is all fixed; "nitrate" is nitrate column when day 0 are raw values
cdata$N<-cdata$n.fixed

#add species names
cdata<-cdata %>%
  mutate(Species=ifelse(species==1, "TT", ifelse(species==2, "CS", ifelse(species==3, "AC", ifelse(species==4, "CH", NA)))))


# plot just the data ------------------ 
cdata %>%
  ggplot(aes(x=day, y=N, color=as.factor(temp))) + geom_point(size=1) + facet_grid(species~temp) +
  geom_errorbar(aes(ymin=N-1.96*(predicted_N_SE), ymax=N+1.96*(predicted_N_SE), 
                    width=.2))
ggsave("N_decay_with_error.pdf")

names(cdata)
# set up models ------------------------


#predict values
rstar_by_temp_pred<-cdata %>% 
  group_by(temp, Species, rep) %>% 
  do(augment(nlsLM(N ~ (a*(day+N_init))/(b+day),
                data= .,  start=list(a=2, b=1),
                control = nls.control(maxiter=100, minFactor=1/204800000))))

cdata %>%
  ggplot(aes(x=day, y=N, group=rep)) + geom_point(size=1) + facet_grid(Species~temp) +
  geom_line(data=rstar_by_temp_pred, aes(y=.fitted), col="red")

#missing a bunch of day 0 values for CS, have a whole bunch of extra 0 values for other species

#fit values
rstar_by_temp<-cdata %>% 
  group_by(temp, Species, rep) %>% 
  do(tidy(nlsLM(N ~ (a*(day+N_init))/(b+day),
                   data= .,  start=list(a=2, b=1),
                   control = nls.control(maxiter=100, minFactor=1/204800000))))

write_csv(rstar_by_temp, "data-processed/rstar_by_temp.csv")

rstar_by_temp<-read_csv("data-processed/rstar_by_temp.csv") %>% 
  filter(term=="a")%>%
  mutate(draw_down=4.4-estimate)

rstar_by_temp %>% 
  #filter(!temp==38 | !Species=="TT")%>%
  #filter(!temp==38 | !Species=="CS")%>%
  #filter(!temp==38 | !Species=="AC")%>%
  ggplot(aes(x = temp, y = estimate, color=as.factor(Species))) + geom_point(size = 2) +
  theme_bw() + facet_grid(~Species) + geom_smooth(method = 'loess') + ylab("Nitrate draw down at equilibrium, uM")
ggsave("figures/rstar_with_temp_2017.pdf", width=7, height=2)


rstar_by_temp_sum<-rstar_by_temp %>% 
  group_by(Species, temp) %>%
  summarize(mean=mean(estimate), error=sd(estimate)/sqrt(n()))

P1<-rstar_by_temp_sum %>% 
  ggplot(aes(x = temp, y = mean, color=as.factor(Species))) + geom_point() +
  theme_bw() + facet_grid(~Species) +
  geom_errorbar(aes(ymin = mean - error, ymax = mean + error), width=0.2) +
  ylab("R-star (uM)")



#read in cell results to display together
cell_results<-read.csv("data-processed/logistic_growth_fits_rep_r-star.csv") 

cell_results_sum<-cell_results %>% 
  mutate(r_corr=r+0.1) %>% #correct for an added death rate of 0.1 per day
  group_by(Species, temp) %>%
  summarize(mean=mean(r_corr), error=sd(r_corr)/sqrt(n()))

P2<-cell_results_sum %>%
  ggplot(aes(x = temp, y = mean, color=as.factor(Species))) + geom_point() +
  theme_bw() + facet_grid(~Species) +
  geom_errorbar(aes(ymin = mean - error, ymax = mean + error), width=0.2) +
  ylab("r (d-1)")

together<-grid.arrange(P1, P2, nrow = 2)
ggsave("figures/2017_direct_estimates.pdf", together, width=8, height=4)


#attempt to show prediction with data

#read in model predictions
best_model_params<-read_csv("data-processed/best_model_params.csv")
best_model_params_booted<-read_csv("data-processed/best_model_params_booted.csv")

minRobs<-rstar_by_temp_sum %>%
  group_by(Species) %>%
  summarize(minR=min(mean))

best_model_params_w_offset<-data.frame()
for (i in unique(best_model_params$Species)){
temp<-best_model_params %>%
  filter(Species==i)  %>%
  mutate(offset=filter(minRobs, Species == i)$minR-min(root))
best_model_params_w_offset<-rbind(best_model_params_w_offset, temp)
}

best_model_params_booted_w_offset<-data.frame()
for (i in unique(best_model_params_booted$Species)){
  temp<-best_model_params_booted %>%
    filter(Species==i)  %>%
    mutate(offset=filter(best_model_params_w_offset, Species == i)$offset[1])
  best_model_params_booted_w_offset<-rbind(best_model_params_booted_w_offset, temp)
}

ggplot() + geom_point(data=rstar_by_temp_sum, aes(x = temp, y = mean, color=as.factor(Species))) +
  theme_bw() + facet_grid(~Species) +  
  geom_errorbar(data=rstar_by_temp_sum, aes(x = temp, ymin = mean - error, ymax = mean + error, color=as.factor(Species)), width=0.2) +
    ylab("R-star (uM)") +
  #geom_line(data=best_model_params_w_offset, aes(x=Temperature, y=root+offset, colour=Species)) + 
  coord_cartesian(ylim = c(1, 5), xlim=c(-1,50)) + 
  geom_ribbon(data=best_model_params_booted_w_offset, 
              aes(x=Temperature, ymin = root_lwr_CI+offset, ymax = root_upr_CI+offset), alpha = .1)

best_model_params_booted_w_offset$offset
