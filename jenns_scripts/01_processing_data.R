# load libraries ----------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)

#### Step 2: create a list of file names for each of the summaries ####

cell_files <- c(list.files("flowcamdata/summaries_040417", full.names = TRUE),
                list.files("flowcamdata/summaries_110417", full.names = TRUE),
                list.files("flowcamdata/summaries_140417", full.names = TRUE),
                list.files("flowcamdata/summaries_180417", full.names = TRUE),
                list.files("flowcamdata/summaries_210417", full.names = TRUE),
                list.files("flowcamdata/summaries_250417", full.names = TRUE),
                list.files("flowcamdata/summaries_280417", full.names = TRUE),
                list.files("flowcamdata/summaries_020517", full.names = TRUE),
                list.files("flowcamdata/summaries_050517", full.names = TRUE),
                list.files("flowcamdata/summaries_090517", full.names = TRUE),
                list.files("flowcamdata/summaries_120517", full.names = TRUE),
                list.files("flowcamdata/summaries_160517", full.names = TRUE),
                list.files("flowcamdata/summaries_190517", full.names = TRUE),
                list.files("flowcamdata/summaries_230517", full.names = TRUE),
                list.files("flowcamdata/summaries_260517", full.names = TRUE),
                list.files("flowcamdata/summaries_300517", full.names = TRUE),
                list.files("flowcamdata/summaries_310517", full.names = TRUE),
                list.files("flowcamdata/summaries_010617", full.names = TRUE),
                list.files("flowcamdata/summaries_020617", full.names = TRUE),
                list.files("flowcamdata/summaries_060617", full.names = TRUE),
                list.files("flowcamdata/summaries_080617", full.names = TRUE),
                list.files("flowcamdata/summaries_080617_highlight", full.names = TRUE),
                list.files("flowcamdata/summaries_090617", full.names = TRUE),
                list.files("flowcamdata/summaries_130617", full.names = TRUE),
                list.files("flowcamdata/summaries_140617", full.names = TRUE),
                list.files("flowcamdata/summaries_150617", full.names = TRUE))
cell_files[1]
names(cell_files) <- cell_files %>% 
  gsub(pattern = ".csv$", replacement = "")

#### Step 3: read in all the files!
#read_csv("flowcamdata/summaries_040417/1-17-2_summary.csv", col_names = FALSE)
        
all_cells <- map_df(cell_files, read_csv, col_names = FALSE, .id = "file_name")

#### Step 4: pull out just the data we want, do some renaming etc.

Rtemp <- all_cells %>% 
  rename(obs_type = X1,
         value = X2) %>% 
  filter(obs_type %in% c("List File", "Start Time", "Particles / ml", "Volume (ABD)")) %>%
  spread(obs_type, value) %>%
  separate(`List File`, into = c("species", "temperature", "replicate"), sep = "-") %>% 
  rename(start_time = `Start Time`,
         cell_density = `Particles / ml`,
         cell_volume = `Volume (ABD)`)


#### Step 5: deal with the date times and separate the replicate names
#set "begin time" as innoculation time
Rtemp$begin.time <- ymd_hms("2017-04-04  11:00:00 AM")
Rtemp$begin.time[Rtemp$species == 2] <- ymd_hms("2017-04-11 11:00:00 AM")
Rtemp$begin.time[Rtemp$species == 3] <- ymd_hms("2017-04-14 11:00:00 AM")
Rtemp$begin.time[Rtemp$species == 4] <- ymd_hms("2017-04-14 11:00:00 AM")

Rtemp$start_time <-  ymd_hms(Rtemp$start_time)

Rtemp$time_since_innoc <- interval(Rtemp$begin.time, Rtemp$start_time)

Rtemp_all <- Rtemp %>% 
  mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
  mutate(time_since_innoc_hours = time_since_innoc/dhours(1)) %>%
  mutate(cell_density = as.numeric(cell_density),
         cell_volume = as.numeric(cell_volume)) %>% 
  mutate(total_biovolume = cell_density * cell_volume) 

names(Rtemp_all)
subset(Rtemp_all, Rtemp_all$species==1 & Rtemp_all$time_since_innoc_days>45)

###fix an error
Rtemp_all$temperature<-ifelse(Rtemp_all$temperature=="07", "10", Rtemp_all$temperature)
#### Step 6: write out the correct csv file!! yay!
write_csv(Rtemp_all, "data-processed/Rtemp_all.csv")
