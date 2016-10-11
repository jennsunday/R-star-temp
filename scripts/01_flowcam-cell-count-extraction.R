library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(lubridate)


#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/Joey/Documents/summaries


#### Step 2: create a list of file names for each of the summaries ####
fnams_aug26 <- list.files("/Users/Joey/Documents/R-star/flowcam-summaries-160826", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
fnams_aug30 <- list.files("/Users/Joey/Documents/R-star/flowcam-summaries-160830", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
fnams_sept6 <- list.files("/Users/Joey/Documents/R-star/flowcam-summaries-160906", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
fnams_sept13 <- list.files("/Users/Joey/Documents/R-star/flowcam-summaries-160913", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file



#### Step 3: create a df with the dataset ID and the cell count ####
Aug_26 <- fnams_aug26 %>% 
	lapply(FUN = function(p) read.csv(p)) %>%
	as.data.frame(.) %>% 
	mutate(List.File = as.character(List.File)) %>% 
dplyr::filter(List.File == "Particles / ml" | List.File == "Start Time") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>% 
	mutate(dataset = rownames(.)) %>%
	mutate(start_time = as.character(V1)) %>%
	mutate(cell_count = as.numeric(as.character(V2))) %>%
	select(-V1) %>% 
	select(-V2)

Aug_30 <- fnams_aug30 %>% 
	lapply(FUN = function(p) read.csv(p)) %>%
	as.data.frame(.) %>%
	mutate(List.File = as.character(List.File)) %>% 
	dplyr::filter(List.File == "Particles / ml" | List.File == "Start Time") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>% 
	mutate(dataset = rownames(.)) %>%
	mutate(start_time = as.character(V1)) %>%
	mutate(cell_count = as.numeric(as.character(V2))) %>%
	select(-V1) %>% 
	select(-V2)



sep6 <- fnams_sept6 %>% 
	lapply(FUN = function(p) read.csv(p)) %>%
	as.data.frame(.) %>%
	mutate(List.File = as.character(List.File)) %>% 
	dplyr::filter(List.File == "Particles / ml" | List.File == "Start Time") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>% 
	mutate(dataset = rownames(.)) %>%
	mutate(start_time = as.character(V1)) %>%
	mutate(cell_count = as.numeric(as.character(V2))) %>%
	select(-V1) %>% 
	select(-V2)

sep13 <- fnams_sept13 %>% 
	lapply(FUN = function(p) read.csv(p)) %>%
	as.data.frame(.) %>%
	mutate(List.File = as.character(List.File)) %>% 
	dplyr::filter(List.File == "Particles / ml" | List.File == "Start Time") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>% 
	mutate(dataset = rownames(.)) %>%
	mutate(start_time = as.character(V1)) %>%
	mutate(cell_count = as.numeric(as.character(V2))) %>%
	select(-V1) %>% 
	select(-V2)

all_points <- bind_rows(Aug_26, sep13, sep6)

all_points_sep <- separate(all_points, dataset, c("species", "replicate"), sep = 2) %>%
	separate(., replicate, c("temperature", "replicate"), sep = 2) %>% 
	separate(., replicate, c("replicate", "file"), sep = 1) %>% 
	select(-file)


#### Step 4: write out the df to a csv ####
write.csv(Rstar_sep_aug26, "data-processed/cell_count_Aug26.csv")


#### Step 5: plot it!

str(all_points_sep)

all_points_sep$start_time <- ymd_hms(all_points_sep$start_time)
all_points_sep %>% 
	filter(species !="CO") %>% 
ggplot(data = ., aes(x = start_time, y = cell_count, group = species, color = species)) + geom_point(size = 2) +
	facet_wrap(~ temperature)
