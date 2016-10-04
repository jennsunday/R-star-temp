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
fnams<- c(list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160816", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160819", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160823", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160826", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160830", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160902", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160906", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160909", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160913", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160916", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160920", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160923", full.names = TRUE), ## find out the names of all the files in data-summary, use full.names to get the relative path for each file
list.files("/Users/Jennifer_Sunday/Documents/R-star/cell_concentrations/flowcam_summaries_160927", full.names = TRUE)) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file


#### Step 3: create a df with the dataset ID and the cell count ####
Augexpt <- fnams %>% 
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


#separate columns into species, temps, and replicates
######
#still need to correct some names and clean this up
######
Rstar_Augexpt <- separate(Augexpt, dataset, c("species", "replicate"), sep = 2) %>% 
	separate(., replicate, c("temperature", "replicate"), sep = 2) %>%
	separate(., replicate, c("replicate", "file"), sep = 1) %>% 
	select(-file)


#make date and time into readable time
str(Rstar_Augexpt)
Rstar_Augexpt$start_time<-ymd_hms(Rstar_Augexpt$start_time)
str(Rstar_sep_aug23)
dim(Rstar_Augexpt)

#factors and numerics
Rstar_Augexpt$species<-as.factor(Rstar_Augexpt$species)
Rstar_Augexpt$temperature<-as.numeric(Rstar_Augexpt$temperature)


#### Step 4: write out the df to a csv ####
write.csv(Rstar_Augexpt, "data-processed/cellcount_Augexpt.csv")


#### Step 5: plot it!
#ggplot(data = Rstar_Augexpt, aes(x = temperature, y = cell_count, group = species, color = species)) + geom_point(size = 4)


par(mfrow=c(1,3))
plot(Rstar_Augexpt$cell_count~Rstar_Augexpt$start_time, ty="n", main="TT")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="TT"& Rstar_Augexpt$temperature==16), points(cell_count~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="TT"& Rstar_Augexpt$temperature==25), points(cell_count~start_time, col=2))
legend("topleft", pch=1, col=c(1,2), c("16°C", "25°C"))

plot(Rstar_Augexpt$cell_count~Rstar_Augexpt$start_time, ty="n", main="CH")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==16), points(cell_count~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="CH"& Rstar_Augexpt$temperature==25), points(cell_count~start_time, col=2))

plot(Rstar_Augexpt$cell_count~Rstar_Augexpt$start_time, ty="n", main="BB")
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==16), points(cell_count~start_time, col=1))
with(subset(Rstar_Augexpt, Rstar_Augexpt$species=="BB"& Rstar_Augexpt$temperature==25), points(cell_count~start_time, col=2))
