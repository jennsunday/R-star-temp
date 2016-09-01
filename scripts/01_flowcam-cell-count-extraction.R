library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)


#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/Joey/Documents/summaries


#### Step 2: create a list of file names for each of the summaries ####
fnams <- list.files("/Users/Joey/Documents/R-star/flowcam-summaries-160826", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file


#### Step 3: create a df with the dataset ID and the cell count ####
Aug_26 <- fnams %>% 
	lapply(FUN = function(p) read.csv(p)) %>%
	as.data.frame(.) %>% 
dplyr::filter(List.File == "Particles / ml") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>% 
	mutate(dataset = rownames(.)) %>%
	mutate(cell_count = as.numeric(as.character(V1))) %>%
	# mutate(biovolume = as.numeric(as.character(V2))) %>% 
	select(-V1) 


Rstar_sep_aug26 <- separate(Aug_26, dataset, c("species", "replicate"), sep = 2) %>%
	separate(., replicate, c("temperature", "replicate"), sep = 2) %>% 
	separate(., replicate, c("replicate", "file"), sep = 1) %>% 
	select(-file)


#### Step 4: write out the df to a csv ####
write.csv(Rstar_sep_aug26, "data-processed/cell_count_Aug26.csv")


#### Step 5: plot it!

ggplot(data = Rstar_sep, aes(x = temperature, y = cell_count, group = species, color = species)) + geom_point()
