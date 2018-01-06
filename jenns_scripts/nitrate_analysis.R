library(tidyr)
library(dplyr)
library(lubridate)
library(plotrix)
library(ggplot2)


nitrate <- read.csv("~/Desktop/data-sheet-nitrate.csv")
nitrate1 <- read.csv("~/Desktop/nitrate-day1.csv")
absorbance <- read.csv("~/Desktop/absorbance_270617.csv")

View(nitrate)

n1 <- nitrate1 %>% 
  select(sample, date.sampled, absorbance.1, absorbance.2) %>% 
  mutate(date = dmy(date.sampled)) %>% 
  filter(!is.na(date)) %>% 
  select(date, sample, absorbance.1) %>% 
  rename(abs = absorbance.1)

n2 <- nitrate %>% 
  select(sample, date.sampled, absorbance_1, absorbance_2) %>% 
  mutate(date = dmy(date.sampled)) %>% 
  filter(!is.na(date)) %>% 
  select(date, sample, absorbance_1) %>% 
  rename(abs = absorbance_1)

n3 <- absorbance %>% 
  select(sample, date.sampled, absorbance.1) %>% 
  mutate(date = dmy(date.sampled)) %>% 
  filter(!is.na(date)) %>% 
  select(date, sample, absorbance.1) %>% 
  rename(abs = absorbance.1)


alln <- bind_rows(n1, n2, n3)

all <- alln %>% 
  filter(!is.na(abs)) %>%
  separate(sample, into = c("species", "temperature", "replicate"))

unique(all$date)



all3 <- all %>% 
  mutate(date2 = as.character(date)) %>% 
  mutate(date3 = ifelse(date2 == "2017-04-13", "2017-04-14", date2)) 


all4 <- all3 %>% 
  mutate(date4 = ifelse(species == "2" & date3 == "2017-04-04", "2017-04-11", date3)) %>% 
  mutate(date = ymd(date4)) %>% 
  select(date, species, temperature, replicate, abs)

write.csv(all4, "~/Desktop/nitrate_data_processed.csv")



all4 %>% 
  group_by(species, temperature, date) %>% 
  summarise_each(funs(mean, std.error), abs) %>% 
  filter(!is.na(mean)) %>% 
  # filter(species == 3) %>% 
  ggplot(aes(x = date, y = mean, color = temperature)) + geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.2) +
  facet_wrap( ~ species) + theme_bw() + ylim(0, 0.1)

all4 %>% 
  unite(uniqueid, species, temperature, replicate, remove = FALSE) %>%
  ggplot(aes(x = date, y = abs, color = temperature, group = uniqueid)) + geom_point() +
  geom_line() +
  facet_wrap( ~ species) + theme_bw() + ylim(0, 0.1)
  
