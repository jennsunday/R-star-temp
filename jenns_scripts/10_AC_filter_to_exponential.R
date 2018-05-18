library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

AC <- read_csv("data/monod_data/AC_15_08_17.csv")

#take raw data, add actual N concentrations
ACN<- AC %>% 	
  mutate(day = Hours.since.Innoc/24) %>% 
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
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


#filter out data where K seems to have been reaACed beyond the variability expected by noise (quite arbitrary)
ACfiltered<- AC %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
  filter(log.Particles.per.ml>2.5) %>%
  filter(  Temperature == 13 & N.Treatment == 8 & day<8 |
             Temperature == 13 & N.Treatment == 7 & day<8|
             Temperature == 13 & N.Treatment == 6 & day<8| 
             Temperature == 13 & N.Treatment == 5 & day<8|
             Temperature == 13 & N.Treatment == 4 & day<8|
             Temperature == 13 & N.Treatment == 3 & day<8 |
             Temperature == 13 & N.Treatment == 2 & day<8|
             Temperature == 13 & N.Treatment == 1 & day<8 |
             Temperature == 13 & N.Treatment == 0 & day<8|
             Temperature == 16 & N.Treatment == 8 & day<8|
             Temperature == 16 & N.Treatment == 7 & day<8|
             Temperature == 16 & N.Treatment == 6 & day<8| 
             Temperature == 16 & N.Treatment == 5 & day<8|
             Temperature == 16 & N.Treatment == 4 & day<8|
             Temperature == 16 & N.Treatment == 3 & day<8|
             Temperature == 16 & N.Treatment == 2 & day<8|
             Temperature == 16 & N.Treatment == 1 & day<8 |
             Temperature == 16 & N.Treatment == 0 & day<8|
             Temperature == 19 & N.Treatment == 8 & day<8|
             Temperature == 19 & N.Treatment == 7 & day<8|
             Temperature == 19 & N.Treatment == 6 & day<8| 
             Temperature == 19 & N.Treatment == 5 & day<8|
             Temperature == 19 & N.Treatment == 4 & day<8 |
             Temperature == 19 & N.Treatment == 3 & day<8 |
             Temperature == 19 & N.Treatment == 2 & day<8 |
             Temperature == 19 & N.Treatment == 1 & day<8 |
             Temperature == 19 & N.Treatment == 0 & day<8|
             Temperature == 22 & N.Treatment == 8 & day<8 |
             Temperature == 22 & N.Treatment == 7 & day<8 |
             Temperature == 22 & N.Treatment == 6 & day<8 | 
             Temperature == 22 & N.Treatment == 5 & day<8 |
             Temperature == 22 & N.Treatment == 4 & day<8 |
             Temperature == 22 & N.Treatment == 3 & day<8 |
             Temperature == 22 & N.Treatment == 2 & day<8 |
             Temperature == 22 & N.Treatment == 1 & day<8 |
             Temperature == 22 & N.Treatment == 0 & day<8|
             Temperature == 25 & N.Treatment == 8 & day<8 |
             Temperature == 25 & N.Treatment == 7 & day<8 |
             Temperature == 25 & N.Treatment == 6 & day<8 | 
             Temperature == 25 & N.Treatment == 5 & day<8 |
             Temperature == 25 & N.Treatment == 4 & day<8 |
             Temperature == 25 & N.Treatment == 3 & day<8 |
             Temperature == 25 & N.Treatment == 2 & day<8 |
             Temperature == 25 & N.Treatment == 1 & day<8 |
             Temperature == 25 & N.Treatment == 0 & day<8|
             Temperature == 28 & N.Treatment == 8 & day<8 |
             Temperature == 28 & N.Treatment == 7 & day<8 |
             Temperature == 28 & N.Treatment == 6 & day<8 | 
             Temperature == 28 & N.Treatment == 5 & day<8 |
             Temperature == 28 & N.Treatment == 4 & day<8 |
             Temperature == 28 & N.Treatment == 3 & day<8 |
             Temperature == 28 & N.Treatment == 2 & day<8 |
             Temperature == 28 & N.Treatment == 1 & day<8 |
             Temperature == 28 & N.Treatment == 0 & day<8)


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

write_csv(ACfilteredN,"data-processed/ACfilteredN.csv")

#Vis linear model over raw data
ACfilteredN %>%
  ggplot(aes(y=log.Particles.per.ml, x=day))  +
  facet_grid(Temperature~N.Treatment) +
  stat_smooth(method=lm) +
  geom_point(data=ACN, color="red", alpha=0.5) + geom_point()
