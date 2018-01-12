library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

CH <- read_csv("data/monod_data/CH_15_08_02.csv")

#take raw data, add actual N concentrations
CHN<- CH %>% 	
  mutate(day = Hours.since.Innoc/24) %>% 
  mutate(N.Treatment = str_replace(N.Treatment, "1", "11")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "2", "22")) %>% 
  mutate(N.Treatment = str_replace(N.Treatment, "3", "33")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "4", "55")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "5", "110")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "6", "220")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "7", "330")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "8", "440")) %>%  
  mutate(N.Treatment = str_replace(N.Treatment, "1105", "55")) %>% 
  mutate(N.Treatment = as.numeric(N.Treatment))

#plot raw data 
CHN %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  #filter(Temperature==16) %>% 
  #mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


#filter out data where K seems to have been reached beyond the variability expected by noise (quite arbitrary)
CHfiltered<- CH %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  filter(Temperature == 13 & N.Treatment == 8 & day>1.5 |
           Temperature == 13 & N.Treatment == 7 & day>1.5 |
           Temperature == 13 & N.Treatment == 6 & day>1.5| 
           Temperature == 13 & N.Treatment == 5 & day>1.5|
           Temperature == 13 & N.Treatment == 4 & day>1.5|
           Temperature == 13 & N.Treatment == 3 & day>0.5|
           Temperature == 13 & N.Treatment == 2 & day>2.5|
           Temperature == 13 & N.Treatment == 1 & day>1.5|
           Temperature == 13 & N.Treatment == 0 & day>1.5|
           Temperature == 16 & N.Treatment == 8 & day<7 & day>0.5|
           Temperature == 16 & N.Treatment == 7 & day<7 & day>0.5|
           Temperature == 16 & N.Treatment == 6 & day<7 & day>0.5| 
           Temperature == 16 & N.Treatment == 5 & day<7 & day>0.5|
           Temperature == 16 & N.Treatment == 4 & day<6 & day>1.5|
           Temperature == 16 & N.Treatment == 3 & day<6 & day>0.5|
           Temperature == 16 & N.Treatment == 2 & day<6 & day>1.5|
           Temperature == 16 & N.Treatment == 1 & day<6 & day>1.5|
           Temperature == 16 & N.Treatment == 0 & day<6 & day>1.5|
           Temperature == 19 & N.Treatment == 8 & day<6|
           Temperature == 19 & N.Treatment == 7 & day<6|
           Temperature == 19 & N.Treatment == 6 & day<6| 
           Temperature == 19 & N.Treatment == 5 & day<6|
           Temperature == 19 & N.Treatment == 4 & day<6|
           Temperature == 19 & N.Treatment == 3 & day<6 |
           Temperature == 19 & N.Treatment == 2 & day<6 |
           Temperature == 19 & N.Treatment == 1 & day<5 & day>1.5|
           Temperature == 19 & N.Treatment == 0 & day<6 & day>1.5|
           Temperature == 22 & N.Treatment == 8 & day<5 |
           Temperature == 22 & N.Treatment == 7 & day<5 |
           Temperature == 22 & N.Treatment == 6 & day<6 | 
           Temperature == 22 & N.Treatment == 5 & day<5 |
           Temperature == 22 & N.Treatment == 4 & day<6 |
           Temperature == 22 & N.Treatment == 3 & day<5 |
           Temperature == 22 & N.Treatment == 2 & day<5 |
           Temperature == 22 & N.Treatment == 1 & day<5 & day>1.5|
           Temperature == 22 & N.Treatment == 0 & day<5 & day>1.5|
           Temperature == 25 & N.Treatment == 8 & day<5 |
           Temperature == 25 & N.Treatment == 7 & day<5 |
           Temperature == 25 & N.Treatment == 6 & day<5 | 
           Temperature == 25 & N.Treatment == 5 & day<5 |
           Temperature == 25 & N.Treatment == 4 & day<4.5 |
           Temperature == 25 & N.Treatment == 3 & day<4 |
           Temperature == 25 & N.Treatment == 2 & day<4 |
           Temperature == 25 & N.Treatment == 1 & day<4 & day>1.5|
           Temperature == 25 & N.Treatment == 0 & day<3.8 & day>1.5|
           Temperature == 28 & N.Treatment == 8 & day<5 |
           Temperature == 28 & N.Treatment == 7 & day<5 |
           Temperature == 28 & N.Treatment == 6 & day<5 | 
           Temperature == 28 & N.Treatment == 5 & day<5 |
           Temperature == 28 & N.Treatment == 4 & day<4 |
           Temperature == 28 & N.Treatment == 3 & day<4 |
           Temperature == 28 & N.Treatment == 2 & day<4 |
           Temperature == 28 & N.Treatment == 1 & day<4 & day>1.5|
           Temperature == 28 & N.Treatment == 0 & day<3.8 & day>1.5)



#take filtered data, add actual N concentrations, plot monod curves
CHfilteredN<- CHfiltered %>% 	
  mutate(N.Treatment = str_replace(N.Treatment, "1", "11")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "2", "22")) %>% 
  mutate(N.Treatment = str_replace(N.Treatment, "3", "33")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "4", "55")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "5", "110")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "6", "220")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "7", "330")) %>%
  mutate(N.Treatment = str_replace(N.Treatment, "8", "440")) %>%  
  mutate(N.Treatment = str_replace(N.Treatment, "1105", "55")) %>% 
  mutate(N.Treatment = as.numeric(N.Treatment))

#plot filtered data - linear model fits...
CHfilteredN %>% 
  #filter(day>0.5) %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+100)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + 
  geom_smooth(method=lm, se=FALSE) 

#plot monod curves with error
linear_r<-CHfilteredN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+100)) %>% 
  filter(day>0.5) %>% 
  group_by(N.Treatment, Temperature) %>% 
  do(tidy(lm(Particles.per.ml ~ day, data=.))) 
linear_r %>% 
  filter(term=="day") %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  theme_bw() + facet_wrap( ~ Temperature) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.2)


write_csv(linear_r, "data-processed/CH_r_lm.csv")
ggsave("figures/CH_monod_curves.png")

#plot linear fit over raw data
linear_r_aug<- CHfilteredN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml+100)) %>% 
  #filter(day>1.5) %>% 
  group_by(N.Treatment, Temperature) %>% 
  do(augment(lm(Particles.per.ml ~ day, data=.))) 

CHN %>% 
  mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  #filter(N.Treatment==33) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + 
  geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + 
  #geom_line(data=subset(linear_r_aug, linear_r_aug$N.Treatment==33), aes(x=day, y=.fitted, colour=factor(N.Treatment)))
  geom_line(data=linear_r_aug, aes(x=day, y=.fitted, colour=factor(N.Treatment)))


ggsave("figures/CH_linear_R_fits.png")


# previous coding working in non-logged data, fitting exponential curve

CHfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ungroup() %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature)
ggsave("figures/CH_monod_curves.png")

CH_r <- CHfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) 

write_csv(CH_r, "data-processed/fitted_r_CH_from_2015.csv")


CHfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ggplot(aes(x = Temperature, y = estimate, color = factor(Temperature))) + geom_point(size = 4) +
  geom_line() + theme_bw() + facet_wrap( ~ N.Treatment)
ggsave("figures/CH_TPC_by_nitrate_curves.png")


#set up the jacknifing
for (j in c(13, 16, 19, 22, 25, 28)){
  for (k in c(0, 11, 22, 33, 55, 110, 220, 440)){
    test<-subset(CHfilteredN, CHfilteredN$Temperature==j & CHfilteredN$N.Treatment==k)
    boot_r<-1:100
    for(i in 1:100){
      mod<-sample_n(test, length(test$day)-1) %>% #take a sample of 1 - total # of days
        do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
                    data= .,  start=list(a=0.01),
                    control = nls.control(maxiter=100, minFactor=1/204800000))))
      boot_r[i]<-mod$estimate
      name<-data.frame(a=unique(boot_r), Temperature=test$Temperature[1], N.Treatment=test$N.Treatment[1])
      assign(paste("CHboot_r_unique", test$Temperature[1], test$N.Treatment[1], sep ="_"), name)
    }
  }
}

# bind all of the objects that start with boot_r_unique
CH_unique_boots<-do.call("rbind", mget(ls(pattern="CHboot_r_unique")))

#plot the jacknifed data - monod
CH_unique_boots %>% 
  group_by(Temperature, N.Treatment) %>% 
  ggplot(aes(x = N.Treatment, y = a, color = factor(N.Treatment))) + geom_point(size = 2) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature)
ggsave("figures/CH_monod_jacknifed.png")