
# load pacakges -----------------------------------------------------------

library(tidyverse)
library(plotrix)
library(stringr)


# read in data ------------------------------------------------------------

aug_data <- 
sept_data <- read_csv("data-processed/cellcount_Septexpt.csv")
oct14 <- read_csv("data-processed/cellcount_Oct14.csv") %>% 
	mutate(temperature = as.integer(temperature))
str(sept_data)
str(oct14)


all_oct <- bind_rows(sept_data, oct14) %>% 
	select(-X1, -X, -unique_ID_1)


data2 <- all_oct %>% 
	separate(start_time, into = c("date", "time"), sep = " ", remove = FALSE) %>% 
	mutate(date_formatted = ymd(date))

sept_data2 <- data2 %>% 
mutate(species = str_replace(species, "ch", "CH")) 


sept_data2$innoc.time <- ymd_hms("2016-09-13 15:19:03")
sept_data2$time_since_innoc <- interval(sept_data2$innoc.time, sept_data2$start_time)
sept_data2 <- sept_data2 %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))


sept_plot <- sept_data2 %>% 
	group_by(species, replicate, temperature) %>% 
	# summarise_each(funs(mean, std.error), cell_count) %>%  
	filter(species != "CO") %>% 
	ggplot(data = ., aes(x = time_since_innoc_hours, y = cell_count, color = species)) + geom_point(size = 2) + 
	# geom_errorbar(aes(ymax = mean + std.error, ymin = mean - std.error)) +
	# geom_line() +
	facet_wrap( ~ temperature) +
	xlab("sample date") +
	ylab("mean cell count") + theme_minimal() + ylim(0, 3*10^5) + xlim(0, 1500)
	
	
	
	
	
	
	group_by(species, date_formatted, temperature) %>% 
	summarise_each(funs(mean, std.error), cell_count) %>%  
	filter(species != "CO") %>% 
	ggplot(data = ., aes(x = date_formatted, y = mean, color = species)) + geom_point(size = 3) + 
	geom_errorbar(aes(ymax = mean + std.error, ymin = mean - std.error)) +
	geom_line() +
	facet_wrap( ~ temperature) +
	xlab("sample date") +
	ylab("mean cell count") + theme_minimal() +
	# scale_y_log10() +
	ylim(0, 3*10^5)

