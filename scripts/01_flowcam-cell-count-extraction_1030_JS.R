library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(lubridate)

##todotoday
#copy new data from this experiment onto my computer
#move summary folders to the summary-only folder
#read in the rest of data

#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/Joey/Documents/summaries
## cp **/*export.csv /Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160913
## cp **/*summary.csv /Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/10_and_30_expt/flowcam_summaries_160913


#### Step 2: create a list of file names for each of the summaries ####
fnams<- c(list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_160913", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_160916", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_160920", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_160923", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_160927", full.names = TRUE),
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_160930", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_161004", full.names = TRUE),
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_161007", full.names = TRUE),
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_161011", full.names = TRUE),
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_161014", full.names = TRUE),
list.files("./data/cell_concentrations/10_and_30/B_flowcam_summaries_161018", full.names = TRUE)) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
## find out the names of all the files in data-summary, use full.names to get the relative path for each file
 ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

#sometimes encounter an error in the following script if some csvs are not correct length, due to csv saving with extra commas (not sure why that occurs)

#### Step 3: create a df with the dataset ID and the cell count ####
Septexpt <- fnams %>% 
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
Rstar_Septexpt <- separate(Septexpt, dataset, c("species", "replicate"), sep = 2) %>% 
	separate(., replicate, c("temperature", "replicate"), sep = 2) %>%
	separate(., replicate, c("replicate", "file"), sep = 1) %>% 
	select(-file)

#### Unite columns to make a unique ID####
Rstar_Septexpt <- Rstar_Septexpt %>%
  unite(unique_ID, species, temperature, replicate, remove = FALSE) 

#make date and time into readable time
str(Rstar_Septexpt)
Rstar_Septexpt$start_time<-ymd_hms(Rstar_Septexpt$start_time)
str(Rstar_Septexpt)


#factors and numerics
Rstar_Septexpt$species<-as.factor(Rstar_Septexpt$species)
Rstar_Septexpt$temperature<-as.numeric(Rstar_Septexpt$temperature)


#### Step 4: write out the df to a csv ####
write.csv(Rstar_Septexpt, "data-processed/cellcount_Septexpt.csv")


#### Step 5: plot it!
#ggplot(data = Rstar_Septexpt, aes(x = temperature, y = cell_count, group = species, color = species)) + geom_point(size = 4)

par(mfrow=c(1,1))
plot(log(Rstar_Septexpt$volume)~Rstar_Septexpt$start_time, type="n", main="TT", xlim=c())
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="TT"& Rstar_Septexpt$temperature==10), points(log(volume)~start_time, col=1))
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="TT"& Rstar_Septexpt$temperature==30), points(log(volume)~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("16°C", "25°C"))

par(mfrow=c(1,3))
plot(log(Rstar_Septexpt$cell_count)~Rstar_Septexpt$start_time, type="n", main="TT", xlim=c())
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="TT"& Rstar_Septexpt$temperature==10), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="TT"& Rstar_Septexpt$temperature==30), points(log(cell_count)~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("16°C", "25°C"))

plot(log(Rstar_Septexpt$cell_count)~Rstar_Septexpt$start_time, ty="n", main="CH")
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="CH"& Rstar_Septexpt$temperature==10), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="CH"& Rstar_Septexpt$temperature==30), points(log(cell_count)~start_time, col=2))

plot(log(Rstar_Septexpt$cell_count)~Rstar_Septexpt$start_time, ty="n", main="BB")
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="BB"& Rstar_Septexpt$temperature==10), points(log(cell_count)~start_time, col=1))
with(subset(Rstar_Septexpt, Rstar_Septexpt$species=="BB"& Rstar_Septexpt$temperature==30), points(log(cell_count)~start_time, col=2))
