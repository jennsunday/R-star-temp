## data cleaning and first plots!
## October 14 2016
## Joey Bernhardt


# load pacakges -----------------------------------------------------------

library(tidyverse)
library(stringr)
library(plotrix)
library(lubridate)
library(gridExtra)



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

aug_data2 <- data1 %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date_formatted = ymd(date))


aug_data2$innoc.time <- ymd_hms("2016-08-12 14:23:57")
aug_data2$time_since_innoc <- interval(aug_data2$innoc.time, aug_data2$start_time)
aug_data2 <- aug_data2 %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))




# plot the mean and std error over time -----------------------------------


aug_plot <- aug_data2 %>% 
	group_by(species, replicate, temperature) %>% 
	# summarise_each(funs(mean, std.error), cell_count) %>%  
	filter(species != "CO") %>% 
	ggplot(data = ., aes(x = time_since_innoc_hours, y = cell_count, color = species)) + geom_point(size = 2) + 
	# geom_errorbar(aes(ymax = mean + std.error, ymin = mean - std.error)) +
	# geom_line() +
				 	facet_wrap( ~ temperature) +
	xlab("sample date") +
	ylab("mean cell count") + theme_minimal() + ylim(0, 3*10^5) + xlim(0, 1500)
	

grid.arrange(aug_plot, sept_plot, ncol = 2)
