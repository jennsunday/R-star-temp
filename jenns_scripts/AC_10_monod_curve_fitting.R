library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

AC <- read_csv("data/monod_data/AC_15_08_17.csv")

#plot raw data 
AC %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  filter(Temperature==22 ) %>% 
  #mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


#filter out data where K seems to have been reaACed beyond the variability expected by noise (quite arbitrary)
ACfiltered<- AC %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  filter(  Temperature == 13 & N.Treatment == 8 & day<7|
           Temperature == 13 & N.Treatment == 7 & day<7|
           Temperature == 13 & N.Treatment == 6 & day<7| 
           Temperature == 13 & N.Treatment == 5 & day<7|
           Temperature == 13 & N.Treatment == 4 & day<7|
           Temperature == 13 & N.Treatment == 3 & day<7|
           Temperature == 13 & N.Treatment == 2 & day<7|
           Temperature == 13 & N.Treatment == 1 & day<7|
           Temperature == 13 & N.Treatment == 0 & day<7|
           Temperature == 16 & N.Treatment == 8 & day<8|
           Temperature == 16 & N.Treatment == 7 & day<8|
           Temperature == 16 & N.Treatment == 6 & day<8| 
           Temperature == 16 & N.Treatment == 5 & day<8|
           Temperature == 16 & N.Treatment == 4 & day<8|
           Temperature == 16 & N.Treatment == 3 & day<8|
           Temperature == 16 & N.Treatment == 2 & day<8|
           Temperature == 16 & N.Treatment == 1 & day<8|
           Temperature == 16 & N.Treatment == 0 & day<8|
           Temperature == 19 & N.Treatment == 8 & day<8|
           Temperature == 19 & N.Treatment == 7 & day<8|
           Temperature == 19 & N.Treatment == 6 & day<8| 
           Temperature == 19 & N.Treatment == 5 & day<8|
           Temperature == 19 & N.Treatment == 4 & day<8|
           Temperature == 19 & N.Treatment == 3 & day<8 |
           Temperature == 19 & N.Treatment == 2 & day<8 |
           Temperature == 19 & N.Treatment == 1 & day<8 |
           Temperature == 19 & N.Treatment == 0 & day<8 |
           Temperature == 22 & N.Treatment == 8 & day<6 |
           Temperature == 22 & N.Treatment == 7 & day<6 |
           Temperature == 22 & N.Treatment == 6 & day<6 | 
           Temperature == 22 & N.Treatment == 5 & day<6 |
           Temperature == 22 & N.Treatment == 4 & day<5.5 |
           Temperature == 22 & N.Treatment == 3 & day<5.5 |
           Temperature == 22 & N.Treatment == 2 & day<5.5 |
           Temperature == 22 & N.Treatment == 1 & day<5.5 |
           Temperature == 22 & N.Treatment == 0 & day<5.5 |
           Temperature == 25 & N.Treatment == 8 & day<6 |
           Temperature == 25 & N.Treatment == 7 & day<6 |
           Temperature == 25 & N.Treatment == 6 & day<6 | 
           Temperature == 25 & N.Treatment == 5 & day<6 |
           Temperature == 25 & N.Treatment == 4 & day<6 |
           Temperature == 25 & N.Treatment == 3 & day<6 |
           Temperature == 25 & N.Treatment == 2 & day<6 |
           Temperature == 25 & N.Treatment == 1 & day<6 |
           Temperature == 25 & N.Treatment == 0 & day<6 |
           Temperature == 28 & N.Treatment == 8 & day<6 |
           Temperature == 28 & N.Treatment == 7 & day<6 |
           Temperature == 28 & N.Treatment == 6 & day<6 | 
           Temperature == 28 & N.Treatment == 5 & day<6 |
           Temperature == 28 & N.Treatment == 4 & day<6 |
           Temperature == 28 & N.Treatment == 3 & day<6 |
           Temperature == 28 & N.Treatment == 2 & day<6 |
           Temperature == 28 & N.Treatment == 1 & day<6 |
           Temperature == 28 & N.Treatment == 0 & day<6  )


#plot filtered data
ACfiltered %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  filter(Temperature==25) %>% 
  #mutate(Particles.per.ml = log(Particles.per.ml)) %>% 
  ggplot(data = ., aes(x = day, y = Particles.per.ml, color = factor(N.Treatment))) + geom_point() +
  facet_wrap( ~ Temperature, scales = "free") + geom_line() + geom_abline(slope = 1, intercept = 5)


#take filtered data, add actual N concentrations, plot monod curves
ACfilteredN<- ACfiltered %>% 	
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


ACfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ungroup() %>% 
  ggplot(aes(x = N.Treatment, y = estimate, color = factor(N.Treatment))) + geom_point(size = 4) +
  geom_line() + theme_bw() + facet_wrap( ~ Temperature) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.1) 
ggsave("figures/AC_monod_curves.png")

AC_r <- ACfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) 

write_ACv(AC_r, "data-processed/fitted_r_AC_from_2015.ACv")


ACfilteredN %>% 
  group_by(Temperature, N.Treatment) %>% 
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
  ggplot(aes(x = Temperature, y = estimate, color = factor(Temperature))) + geom_point(size = 1) +
  geom_line() + theme_bw() + facet_wrap( ~ N.Treatment) + 
  geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate+std.error), width=.1) 
ggsave("figures/AC_TPC_by_nitrate_curves.png")


#set up the bootstrapping
for (j in c(13, 16, 19, 22, 25, 28)){
for (k in c(0, 11, 22, 33, 55, 110, 440)){
test<-subset(ACfilteredN, ACfilteredN$Temperature==13 & ACfilteredN$N.Treatment==k)
boot_r<-1:100
for(i in 1:100){
  mod<-sample_n(test, length(test$day)-1) %>%
  do(tidy(nls(Particles.per.ml ~ 75 * (1+a)^(Hours.since.Innoc),
              data= .,  start=list(a=0.01),
              control = nls.control(maxiter=100, minFactor=1/204800000))))
  boot_r[i]<-mod$estimate
  name<-data.frame(a=unique(boot_r), Temperature=test$Temperature[1], N.Treatment=test$N.Treatment[1])
  assign(paste("boot_r_unique", test$Temperature[1], test$N.Treatment[1], sep ="_"), name)
 }
}
}
bootedr<-rbind(boot_r_unique_13_0, 
               boot_r_unique_13_11, 
               boot_r_unique_13_22, 
             
               , 
               boot_r_unique_13_11, 
               boot_r_unique_13_22,boot_r_unique_22, boot_r_unique_25, boot_r_unique_28)

grep("boot_r_unique", R_GlobalEnv)

rbind("boot_r_uniqu_e*")
boot_r_unique_13_11

list(boot_r_unique*)
bootedr %>% 
  group_by(Temperature) %>% 
  ggplot(aes(x = Temperature, y = r, color = factor(Temperature))) + geom_point(size = 1) +
  geom_line() + theme_bw() + facet_wrap( ~ N.Treatment) 