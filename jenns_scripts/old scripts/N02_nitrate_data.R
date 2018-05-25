
#Goal bring together raw nitrate absorbance data from multiple files

data<-rbind(read.csv("nitrate_data/absorbance_data_raw/nitrate_data_Apr04.csv"),
            read.csv("nitrate_data/absorbance_data_raw/nitrate_data_Apr11.csv"),
            read.csv("nitrate_data/absorbance_data_raw/nitrate_data_Apr14.csv"),
            read.csv("nitrate_data/absorbance_data_raw/nitrate_data_Apr18.csv"), 
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_Apr25.csv"),
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_May02.csv"),
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_May09.csv"),
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_May12.csv"),
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_May16.csv"),
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_May19.csv"),
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_May23.csv"),
          read.csv("nitrate_data/absorbance_data_raw/nitrate_data_enddays.csv"))
dim(data)
head(data) 
tail(data)

#write csv
write.csv(data, "data-processed/nitratedatacombined.csv")

