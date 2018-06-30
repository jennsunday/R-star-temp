library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)
library(dplyr)

CS <- read_csv("data/monod_data/CS_15_08_04.csv")

CS<-CS %>% 
  select(-Hours.since.Innoc, -corrected_hour) %>% 
  rename(Hours.since.Innoc=corrected.hours.since.innoc)


#take raw data, add actual N concentrations
CSN<- CS %>% 	
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



#filter out data where K seems to have been reaCSed beyond the variability expected by noise (quite arbitrary)
CSfiltered<- CS %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  #filter(day>0.5) %>% 
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
  filter(  Temperature == 13 & N.Treatment == 8 & day<8 & day>0.5|
             Temperature == 13 & N.Treatment == 7 & day<8 & day>0.5|
             Temperature == 13 & N.Treatment == 6 & day<8 & day>0.5| 
             Temperature == 13 & N.Treatment == 5 & day<8 & day>0.5|
             Temperature == 13 & N.Treatment == 4 & day<8 & day>0.5|
             Temperature == 13 & N.Treatment == 3 & day<8 & day>0.5|
             Temperature == 13 & N.Treatment == 2 & day<8 & day>0.5|
             Temperature == 13 & N.Treatment == 1 & day<8 & day>0.5|
             Temperature == 13 & N.Treatment == 0 & day<8|
             Temperature == 16 & N.Treatment == 8 & day<6|
             Temperature == 16 & N.Treatment == 7 & day<6|
             Temperature == 16 & N.Treatment == 6 & day<6| 
             Temperature == 16 & N.Treatment == 5 & day<6|
             Temperature == 16 & N.Treatment == 4 & day<6|
             Temperature == 16 & N.Treatment == 3 & day<6|
             Temperature == 16 & N.Treatment == 2 & day<6|
             Temperature == 16 & N.Treatment == 1 & day<6|
             Temperature == 16 & N.Treatment == 0 & day<8|
             Temperature == 19 & N.Treatment == 8 & day<6|
             Temperature == 19 & N.Treatment == 7 & day<6|
             Temperature == 19 & N.Treatment == 6 & day<6| 
             Temperature == 19 & N.Treatment == 5 & day<6|
             Temperature == 19 & N.Treatment == 4 & day<6|
             Temperature == 19 & N.Treatment == 3 & day<6 |
             Temperature == 19 & N.Treatment == 2 & day<6 |
             Temperature == 19 & N.Treatment == 1 & day<6 |
             Temperature == 19 & N.Treatment == 0 & day<8|
             Temperature == 22 & N.Treatment == 8 & day<6 |
             Temperature == 22 & N.Treatment == 7 & day<6 |
             Temperature == 22 & N.Treatment == 6 & day<6 | 
             Temperature == 22 & N.Treatment == 5 & day<6 |
             Temperature == 22 & N.Treatment == 4 & day<5.5 |
             Temperature == 22 & N.Treatment == 3 & day<5.5 |
             Temperature == 22 & N.Treatment == 2 & day<5.5 |
             Temperature == 22 & N.Treatment == 1 & day<5.5 |
             Temperature == 22 & N.Treatment == 0 & day<8|
             Temperature == 25 & N.Treatment == 8 & day<6 |
             Temperature == 25 & N.Treatment == 7 & day<6 |
             Temperature == 25 & N.Treatment == 6 & day<6 | 
             Temperature == 25 & N.Treatment == 5 & day<6 |
             Temperature == 25 & N.Treatment == 4 & day<6 |
             Temperature == 25 & N.Treatment == 3 & day<5.5 |
             Temperature == 25 & N.Treatment == 2 & day<5.5 |
             Temperature == 25 & N.Treatment == 1 & day<4.5 |
             Temperature == 25 & N.Treatment == 0 & day<8|
             Temperature == 28 & N.Treatment == 8 & day<6 |
             Temperature == 28 & N.Treatment == 7 & day<6 |
             Temperature == 28 & N.Treatment == 6 & day<6 | 
             Temperature == 28 & N.Treatment == 5 & day<6 |
             Temperature == 28 & N.Treatment == 4 & day<6 |
             Temperature == 28 & N.Treatment == 3 & day<6 |
             Temperature == 28 & N.Treatment == 2 & day<6 |
             Temperature == 28 & N.Treatment == 1 & day<5.5 |
             Temperature == 28 & N.Treatment == 0 & day<8)

#take filtered data, add actual N concentrations
CSfilteredN <- CSfiltered %>% 	
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


write_csv(CSfilteredN,"data-processed/CSfilteredN.csv")
write_csv(CSN,"data-processed/CSN.csv")

CSN<-read_csv("data-processed/CSN.csv")
CSfilteredN<-read_csv("data-processed/CSfilteredN.csv")

#Vis linear model over raw data
library(cowplot)
CSfilteredN %>%
  ggplot(aes(y=log.Particles.per.ml, x=day))  +
  facet_grid(Temperature~N.Treatment) +
  stat_smooth(method=lm) +
  geom_point(data=TTN, aes(y=log.Particles.per.ml), color=grey(0.3), alpha=0.5) +
  geom_point() + labs(y="Log cell density", x="time, d") + 
  ggtitle("CS")
ggsave("figures/CS_raw_data.pdf", width=8, height=7)  
