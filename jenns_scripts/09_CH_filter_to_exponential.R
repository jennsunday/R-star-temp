library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

CH <- read_csv("data/monod_data/CH_15_08_02.csv")

names(CH)
#take raw data, add actual N concentrations
CHN<- CH %>% 	
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

#
CHfiltered<- CH %>% 
  mutate(day = Hours.since.Innoc/24) %>% 
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	
  filter(Temperature == 13 & N.Treatment == 8 |
           Temperature == 13 & N.Treatment == 7 |
           Temperature == 13 & N.Treatment == 6 | 
           Temperature == 13 & N.Treatment == 5 |
           Temperature == 13 & N.Treatment == 4 |
           Temperature == 13 & N.Treatment == 3 |
           Temperature == 13 & N.Treatment == 2 |
           Temperature == 13 & N.Treatment == 1 |
           Temperature == 13 & N.Treatment == 0 |
           Temperature == 16 & N.Treatment == 8 & day<7 |
           Temperature == 16 & N.Treatment == 7 & day<7 |
           Temperature == 16 & N.Treatment == 6 & day<7 | 
           Temperature == 16 & N.Treatment == 5 & day<7 |
           Temperature == 16 & N.Treatment == 4 & day<6 |
           Temperature == 16 & N.Treatment == 3 & day<6 |
           Temperature == 16 & N.Treatment == 2 & day<6 |
           Temperature == 16 & N.Treatment == 1 & day<6 |
           Temperature == 16 & N.Treatment == 0 & day<8 |
           Temperature == 19 & N.Treatment == 8 & day<6|
           Temperature == 19 & N.Treatment == 7 & day<6|
           Temperature == 19 & N.Treatment == 6 & day<6| 
           Temperature == 19 & N.Treatment == 5 & day<6|
           Temperature == 19 & N.Treatment == 4 & day<6|
           Temperature == 19 & N.Treatment == 3 & day<6 |
           Temperature == 19 & N.Treatment == 2 & day<5 |
           Temperature == 19 & N.Treatment == 1 & day<5 |
           Temperature == 19 & N.Treatment == 0 & day<8 |
           Temperature == 22 & N.Treatment == 8 & day<6 |
           Temperature == 22 & N.Treatment == 7 & day<6 |
           Temperature == 22 & N.Treatment == 6 & day<6 | 
           Temperature == 22 & N.Treatment == 5 & day<6 |
           Temperature == 22 & N.Treatment == 4 & day<6 |
           Temperature == 22 & N.Treatment == 3 & day<5 |
           Temperature == 22 & N.Treatment == 2 & day<5 |
           Temperature == 22 & N.Treatment == 1 & day<3.5|
           Temperature == 22 & N.Treatment == 0 & day<3.5|
           Temperature == 25 & N.Treatment == 8 & day<4.5 |
           Temperature == 25 & N.Treatment == 7 & day<4.5 |
           Temperature == 25 & N.Treatment == 6 & day<4.5 | 
           Temperature == 25 & N.Treatment == 5 & day<4.5 |
           Temperature == 25 & N.Treatment == 4 & day<4 |
           Temperature == 25 & N.Treatment == 3 & day<4 |
           Temperature == 25 & N.Treatment == 2 & day<4 & day>0.8|
           Temperature == 25 & N.Treatment == 1 & day<3.5|
           Temperature == 25 & N.Treatment == 0 & day<3.5|
           Temperature == 28 & N.Treatment == 8 & day<4.5 |
           Temperature == 28 & N.Treatment == 7 & day<4.5 |
           Temperature == 28 & N.Treatment == 6 & day<4.5 | 
           Temperature == 28 & N.Treatment == 5 & day<4 |
           Temperature == 28 & N.Treatment == 4 & day<4 |
           Temperature == 28 & N.Treatment == 3 & day<4 |
           Temperature == 28 & N.Treatment == 2 & day<4 |
           Temperature == 28 & N.Treatment == 1 & day<3.5|
           Temperature == 28 & N.Treatment == 0 & day<3.5)



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

write_csv(CHfilteredN,"data-processed/CHfilteredN.csv")
write_csv(CHN,"data-processed/CHN.csv")

CHN<-read_csv("data-processed/CHN.csv")
CHfilteredN<-read_csv("data-processed/CHfilteredN.csv")

#Vis linear model over raw data
library(cowplot)
CHfilteredN %>%
  ggplot(aes(y=log.Particles.per.ml, x=day))  +
  facet_grid(Temperature~N.Treatment) +
  stat_smooth(method=lm) +
  geom_point(data=TTN, aes(y=log.Particles.per.ml), color=grey(0.3), alpha=0.5) +
  geom_point() + labs(y="Log cell density", x="time, d") + 
  ggtitle("CH")
ggsave("figures/CH_raw_data.pdf", width=8, height=7)  
