#Goal: create functions to fit logistic growth curve to nitrate data

#read in libraries
library(simecol)
library(tidyverse)
library(broom)
library(minpack.lm)

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



rstar_by_temp %>% 
  filter(term=="a") %>% 
  filter(!temp==38 | !Species=="TT")%>%
  filter(!temp==38 | !Species=="CS")%>%
  filter(!temp==38 | !Species=="AC")%>%
  ggplot(aes(x = temp, y = estimate, color=as.factor(Species))) + geom_point(size = 2) +
  theme_bw() + facet_grid(~Species) + geom_smooth(method = 'loess')
ggsave("rstar_with_temp_2017.pdf", width=7, height=2)
