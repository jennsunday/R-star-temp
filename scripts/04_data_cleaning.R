## data cleaning and first plots!
## October 14 2016
## Joey Bernhardt


# load pacakges -----------------------------------------------------------

library(tidyverse)
library(stringr)
library(plotrix)
library(lubridate)



# read in data ------------------------------------------------------------

data <- read_csv("data-processed/cellcount_Augexpt.csv") %>% 
	select(-X1)


# clean up the weirdo species names ---------------------------------------

data1 <- data %>% 
	mutate(temperature = str_replace(temperature, "26|53|52", "25")) %>% 
	mutate(species = str_replace(species, "B2", "BB")) %>% 
	mutate(species = str_replace(species, "BL", "CO")) %>% 
	mutate(species = str_replace(species, "ch", "CH")) %>% 
	mutate(species = str_replace(species, "X2", "CO"))
	

unique(data$species)


# separate the date-time column into a date only col ----------------------

data2 <- data1 %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date_formatted = ymd(date))



# plot the mean and std error over time -----------------------------------


data2 %>% 
	group_by(species, date_formatted, temperature) %>% 
	summarise_each(funs(mean, std.error), cell_count) %>%  
	filter(species != "CO") %>% 
	ggplot(data = ., aes(x = date_formatted, y = mean, color = species)) + geom_point(size = 2) + 
	geom_errorbar(aes(ymax = mean + std.error, ymin = mean - std.error)) +
	geom_line() +
				 	facet_wrap( ~ temperature) +
	xlab("sample date") +
	ylab("mean cell count") + theme_minimal()
	
