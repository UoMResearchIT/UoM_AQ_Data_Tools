set.seed(171219)
#-----------------------------------------------------------------------------------
# Author:       Manuele Reani
# Date:         04/03/2020
# Institution:  SME - CUHK, Shenzhen (China)
# Object:       Pollen and Weather Data Cleaning
#-----------------------------------------------------------------------------------


# set working dir
setwd("/mnt/iusers01/support/mbessdl2/scratch/Turing_AQ_Data_Prep")

original_filename <- "midas_source.csv"
output_filename <- "station_data_clean2.csv"

# -------------------------------------------------------------------
# stations data
# -------------------------------------------------------------------
message("reading original station data")
station_data <- read.csv(file=original_filename,header=TRUE,sep=",")

message("stripping out unneeded columns, and naming required columns")
# only keep relevant features for station
station_data <- station_data[,c(1,3,4,8)]
# name these features 
names(station_data) <- 
  c("Station", "Latitude", "Longitude", "Postcode")
  
message("dropping empty records")
# eliminate empty record
station_data <- station_data[-c(29870,26186,26190,26185, 26189),]

message("stripping out the primary postcode area")
# keep only the postcode area
station_data$Postcode <- gsub("[^a-zA-Z]", "", station_data$Postcode)

message("writing required station records to output file")
# save this data ############################################################
write.csv(station_data, file = output_filename, row.names=F)
