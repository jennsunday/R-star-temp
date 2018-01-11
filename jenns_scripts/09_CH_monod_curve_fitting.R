library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

CH <- read_csv("data/monod_data/CH_15_08_02.csv")

#plot raw data 
CH %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  filter(Temperature==16) %>% 
  #mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


#filter out data where K seems to have been reached beyond the variability expected by noise (quite arbitrary)
CHfiltered<- CH %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  filter(Temperature == 13 & N.Treatment == 8 |
           Temperature == 13 & N.Treatment == 7 |
           Temperature == 13 & N.Treatment == 6 | 
           Temperature == 13 & N.Treatment == 5 |
           Temperature == 13 & N.Treatment == 4 |
           Temperature == 13 & N.Treatment == 3 |
           Temperature == 13 & N.Treatment == 2 |
           Temperature == 13 & N.Treatment == 1 |
           Temperature == 13 & N.Treatment == 0 |
           Temperature == 16 & N.Treatment == 8 & day<6|
           Temperature == 16 & N.Treatment == 7 & day<6|
           Temperature == 16 & N.Treatment == 6 & day<6| 
           Temperature == 16 & N.Treatment == 5 & day<6|
           Temperature == 16 & N.Treatment == 4 & day<6 |
           Temperature == 16 & N.Treatment == 3 & day<6 |
           Temperature == 16 & N.Treatment == 2 & day<6 |
           Temperature == 16 & N.Treatment == 1 & day<6 |
           Temperature == 16 & N.Treatment == 0 & day<6 |
           Temperature == 19 & N.Treatment == 8 & day<6|
           Temperature == 19 & N.Treatment == 7 & day<6|
           Temperature == 19 & N.Treatment == 6 & day<6| 
           Temperature == 19 & N.Treatment == 5 & day<6|
           Temperature == 19 & N.Treatment == 4 & day<6|
           Temperature == 19 & N.Treatment == 3 & day<6 |
           Temperature == 19 & N.Treatment == 2 & day<6 |
           Temperature == 19 & N.Treatment == 1 & day<6 |
           Temperature == 19 & N.Treatment == 0 & day<6 |
           Temperature == 22 & N.Treatment == 8 & day<5 |
           Temperature == 22 & N.Treatment == 7 & day<5 |
           Temperature == 22 & N.Treatment == 6 & day<5 | 
           Temperature == 22 & N.Treatment == 5 & day<5 |
           Temperature == 22 & N.Treatment == 4 & day<5 |
           Temperature == 22 & N.Treatment == 3 & day<5 |
           Temperature == 22 & N.Treatment == 2 & day<5 |
           Temperature == 22 & N.Treatment == 1 & day<5 |
           Temperature == 22 & N.Treatment == 0 & day<5 |
           Temperature == 25 & N.Treatment == 8 & day<4 |
           Temperature == 25 & N.Treatment == 7 & day<4 |
           Temperature == 25 & N.Treatment == 6 & day<4 | 
           Temperature == 25 & N.Treatment == 5 & day<4 |
           Temperature == 25 & N.Treatment == 4 & day<4 |
           Temperature == 25 & N.Treatment == 3 & day<4 |
           Temperature == 25 & N.Treatment == 2 & day<4 |
           Temperature == 25 & N.Treatment == 1 & day<4 |
           Temperature == 25 & N.Treatment == 0 & day<3.8 |
           Temperature == 28 & N.Treatment == 8 & day<4 |
           Temperature == 28 & N.Treatment == 7 & day<4 |
           Temperature == 28 & N.Treatment == 6 & day<4 | 
           Temperature == 28 & N.Treatment == 5 & day<4 |
           Temperature == 28 & N.Treatment == 4 & day<4 |
           Temperature == 28 & N.Treatment == 3 & day<4 |
           Temperature == 28 & N.Treatment == 2 & day<4 |
           Temperature == 28 & N.Treatment == 1 & day<4 |
           Temperature == 28 & N.Treatment == 0 & day<3.8  )


#plot filtered data
CHfiltered %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  #mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


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
