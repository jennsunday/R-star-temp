library(tidyverse)
library(broom)
library(minpack.lm)
library(purrr)
library(stringr)

TT <- read_csv("data/monod_data/TT_15_08_24.csv")


#take raw data, add actual N concentrations
TTN<- TT %>% 	
  filter(N.Treatment!=5 | Temperature!=28 | MINUTE!=22) %>% 
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

#filter out data where K seems to have been reached beyond the variability expected by noise (quite arbitrary)
TTfiltered<- TT %>% 
  mutate(day = Hours.since.Innoc/24) %>% #convert hours to day
  filter(N.Treatment!=5 | Temperature!=28 | MINUTE!=22) %>%  #remove (why?)
  mutate(log.Particles.per.ml = log(Particles.per.ml+1)) %>% 	#create a column with log(cell density)
  filter( Temperature == 13 & N.Treatment == 8 |
            Temperature == 13 & N.Treatment == 7 |
            Temperature == 13 & N.Treatment == 6 | 
            Temperature == 13 & N.Treatment == 5 |
            Temperature == 13 & N.Treatment == 4 |
            Temperature == 13 & N.Treatment == 3 |
            Temperature == 13 & N.Treatment == 2 |
            Temperature == 13 & N.Treatment == 1 & day<6 |
            Temperature == 13 & N.Treatment == 0  |
            Temperature == 16 & N.Treatment == 8 |
            Temperature == 16 & N.Treatment == 7 |
            Temperature == 16 & N.Treatment == 6 | 
            Temperature == 16 & N.Treatment == 5  |
            Temperature == 16 & N.Treatment == 4 |
            Temperature == 16 & N.Treatment == 3 |
            Temperature == 16 & N.Treatment == 2 |
            Temperature == 16 & N.Treatment == 1 & day<6 |
            Temperature == 16 & N.Treatment == 0 |
            Temperature == 19 & N.Treatment == 8  |
            Temperature == 19 & N.Treatment == 7  |
            Temperature == 19 & N.Treatment == 6 | 
            Temperature == 19 & N.Treatment == 5 |
            Temperature == 19 & N.Treatment == 4 |
            Temperature == 19 & N.Treatment == 3 |
            Temperature == 19 & N.Treatment == 2 |
            Temperature == 19 & N.Treatment == 1 |
            Temperature == 19 & N.Treatment == 0 |
            Temperature == 22 & N.Treatment == 8 | 
            Temperature == 22 & N.Treatment == 7 |
            Temperature == 22 & N.Treatment == 6 | 
            Temperature == 22 & N.Treatment == 5 |
            Temperature == 22 & N.Treatment == 4 |
            Temperature == 22 & N.Treatment == 3 & day<4.5|
            Temperature == 22 & N.Treatment == 2 |
            Temperature == 22 & N.Treatment == 1 |
            Temperature == 22 & N.Treatment == 0 | 		
            Temperature == 25 & N.Treatment == 8 |
            Temperature == 25 & N.Treatment == 7 |
            Temperature == 25 & N.Treatment == 6 | 
            Temperature == 25 & N.Treatment == 5 |
            Temperature == 25 & N.Treatment == 4 & day<3.5 |
            Temperature == 25 & N.Treatment == 3 & day<3.5 |
            Temperature == 25 & N.Treatment == 2 |
            Temperature == 25 & N.Treatment == 1 & day<3.5|
            Temperature == 25 & N.Treatment == 0 | 		
            Temperature == 28 & N.Treatment == 8 |
            Temperature == 28 & N.Treatment == 7 |
            Temperature == 28 & N.Treatment == 6 | 
            Temperature == 28 & N.Treatment == 5 |
            Temperature == 28 & N.Treatment == 4 & day<3 |
            Temperature == 28 & N.Treatment == 3 & day<3 |
            Temperature == 28 & N.Treatment == 2 |
            Temperature == 28 & N.Treatment == 1 & day<4|
            Temperature == 28 & N.Treatment == 0 ) 


#take filtered data, add actual N concentrations, plot monod curves
TTfilteredN<- TTfiltered %>% 	
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

write_csv(TTfilteredN,"data-processed/TTfilteredN.csv")
write_csv(TTN,"data-processed/TTN.csv")
TTN<-read_csv("data-processed/TTN.csv")
TTfilteredN<-read_csv("data-processed/TTfilteredN.csv")


#Vis linear model over raw data
library(cowplot)

TTfilteredN %>%
  ggplot(aes(y=log.Particles.per.ml, x=day))  +
  facet_grid(Temperature~N.Treatment) +
  stat_smooth(method=lm) +
  geom_point(data=TTN, aes(y=log.Particles.per.ml), color=grey(0.3), alpha=0.5) +
  geom_point() + labs(y="Log cell density", x="time, d") + 
  ggtitle("TT")
ggsave("figures/TT_raw_data.pdf", width=8, height=7)  
