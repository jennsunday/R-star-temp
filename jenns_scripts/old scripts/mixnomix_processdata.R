# load libraries ----------------------------------------------------------

library(tidyverse)
library(purrr)
library(lubridate)
library(stringr)

#### Step 2: create a list of file names for each of the summaries ####

mixcell_files <- c(list.files("flowcamdata/summaries_mix_310517", full.names = TRUE),
                   list.files("flowcamdata/summaries_mix_010617", full.names = TRUE))

mixcell_files[1]
names(mixcell_files) <- mixcell_files %>% 
  gsub(pattern = ".csv$", replacement = "")

#### Step 3: read in all the files!
#read_csv("flowcamdata/summaries_040417/1-17-2_summary.csv", col_names = FALSE)

mixall_cells <- map_df(mixcell_files, read_csv, col_names = FALSE, .id = "file_name")

#### Step 4: pull out just the data we want, do some renaming etc.

mixnomix <- mixall_cells %>% 
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
mixnomix$begin.time <- ymd_hms("2017-04-04  11:00:00 AM")
mixnomix$begin.time[mixnomix$species == "mx2"] <- ymd_hms("2017-04-11 11:00:00 AM")
mixnomix$begin.time[mixnomix$species == "mx3"] <- ymd_hms("2017-04-14 11:00:00 AM")
mixnomix$begin.time[mixnomix$species == "mx4"] <- ymd_hms("2017-04-14 11:00:00 AM")
mixnomix$begin.time[mixnomix$species == "nm1"] <- ymd_hms("2017-04-04  11:00:00 AM")
mixnomix$begin.time[mixnomix$species == "nm2"] <- ymd_hms("2017-04-11 11:00:00 AM")
mixnomix$begin.time[mixnomix$species == "nm3"] <- ymd_hms("2017-04-14 11:00:00 AM")
mixnomix$begin.time[mixnomix$species == "nm4"] <- ymd_hms("2017-04-14 11:00:00 AM")

mixnomix$start_time <-  ymd_hms(mixnomix$start_time)

mixnomix$time_since_innoc <- interval(mixnomix$begin.time, mixnomix$start_time)

mixnomix_all <- mixnomix %>% 
  mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
  mutate(time_since_innoc_hours = time_since_innoc/dhours(1)) %>%
  mutate(cell_density = as.numeric(cell_density),
         cell_volume = as.numeric(cell_volume)) %>% 
  mutate(total_biovolume = cell_density * cell_volume) 

names(mixnomix_all)

#### Step 6: write out the correct csv file!! yay!
write_csv(mixnomix_all, "data-processed/mixnomix_all.csv")
