#

# load packages -----------------------------------------------------------

library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(lubridate)

#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/Joey/Documents/summaries
## cp **/*export.csv /Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160913
## cp **/*summary.csv /Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/10_and_30_expt/flowcam_summaries_160913


#### Step 2: create a list of file names for each of the summaries ####
fnams <- list.files("/Users/Joey/Documents/R-star/flowcam-summaries-161014", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
## find out the names of all the files in data-summary, use full.names to get the relative path for each file
fnams

#### Step 3: create a df with the dataset ID and the cell count ####
Oct_14 <- fnams %>% 
	lapply(FUN = function(p) read.csv(p)) %>% 
	as.data.frame(.) %>% 
	mutate(List.File = as.character(List.File)) %>% 
	dplyr::filter(List.File == "Particles / ml" | List.File == "Start Time" | List.File == "Volume (ABD)") %>% 
	select(- starts_with("List")) %>%
	t(.) %>% 
	as.data.frame() %>% 
	mutate(dataset = rownames(.)) %>% 
	mutate(start_time = as.character(V1)) %>%
	mutate(cell_count = as.numeric(as.character(V2))) %>% 
	mutate(volume = as.numeric(as.character(V3))) %>% 
	select(-V1) %>% 
	select(-V2) %>% 
	select(-V3)

######
#still need to correct some names and clean this up
######

#separate columns into species, temps, and replicates
Rstar_Oct14 <- separate(Oct_14, dataset, c("species", "replicate"), sep = 2) %>% 
	separate(., replicate, c("temperature", "replicate"), sep = 2) %>%
	separate(., replicate, c("replicate", "file"), sep = 1) %>% 
	select(-file)

#### Unite columns to make a unique ID####
Rstar_Oct14 <- Rstar_Oct14 %>%
	unite(unique_ID, species, temperature, replicate, remove = FALSE) 

#make date and time into readable time

Rstar_Oct14$start_time<-ymd_hms(Rstar_Oct14$start_time)



# #factors and numerics
# Rstar_Septexpt$species<-as.factor(Rstar_Septexpt$species)
# Rstar_Septexpt$temperature<-as.numeric(Rstar_Septexpt$temperature)


#### Step 4: write out the df to a csv ####
write.csv(Rstar_Oct14, "data-processed/cellcount_Oct14.csv")
