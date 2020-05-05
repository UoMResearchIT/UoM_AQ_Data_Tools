set.seed(171219)
#-----------------------------------------------------------------------------------
# Author:       Manuele Reani
# Date:         04/03/2020
# Institution:  SME - CUHK, Shenzhen (China)
# Object:       Pollen and Weather Data Cleaning
#-----------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# FUNCTION:     loadPackages(package.args)
# INPUT:        vector
# OUTPUT:       void
# DESCRIPTION:  Loads required packages.
#                
#---------------------------------------------------------------------------------
loadPackages <- function(package.args)
{ 
  for(i in package.args)
  {
    if(!is.element(i, .packages(all.available = TRUE)))
    {
      cat("\nPackage <", i, "> not found, attempting to add it...")
      install.packages(i)
    }
    library(i, character.only = TRUE)
  }
}
#---------------------------------------------------------------------------------
# FUNCTION:     initialize()
# INPUT:        void
# OUTPUT:       void
# DESCRIPTION:  Set up function for adding packages and other source data
#               
#---------------------------------------------------------------------------------
initialize <- function()
{
  # load packages
#  package.args <- c("car","msm", "png","jpeg", "gplots", "dplyr", "plyr","tidyr", 
#                    "matrixStats",  "lattice","ggplot2", "gtools", 
#                    "dbscan", "stringdist", "utils", "qualV", "stringi", "dplyr", 
#                    "stringr", "rjson","lsmeans","multcomp","lme4","nlme","MuMIn",
#                    "effsize","heplots","DescTools","irr","reshape","psych","effsize",
#                    "ggthemes","ggmap","maps","mapdata","mice","VIM")
  package.args <- c("dplyr")
  loadPackages(package.args)
}

#initialize()

message("finished initialization")

used_MICE <- FALSE

station_dataset_in <- "station_data_clean.csv"


# set working dir
setwd("/mnt/iusers01/support/mbessdl2/scratch/Turing_AQ_Data_Prep")

#----------------------------------------------------------------------------
#  Choose input files based on if MICE was used or not.
#----------------------------------------------------------------------------
if (used_MICE){
	message("Using input data with MICE imputation")
	# input file names
	weather_dataset_in <- "weather_df_Final_comp.csv"
	pollen_dataset_in <- "pollen_df_Final_comp.csv"
	rain_dataset_in <- "rain_df_Final_comp.csv"
	# output file names
	weather_output_file <- "out_weather.csv"
	pollen_output_file <- "out_pollen.csv"
	rain_output_file <- "out_rain.csv"
} else {
	message("Using plain input data")
	# input file names
	weather_dataset_in <- "weather_df_Max_Mean_Intermediate.csv"
	pollen_dataset_in <- "pollen_df_Max_Mean_Intermediate.csv"
	rain_dataset_in <- "rain_df_Max_Mean_Intermediate.csv"
	# output file names
	weather_output_file <- "Raw_Weather_Data.csv"
	pollen_output_file <- "Raw_Pollen_Data.csv"
	rain_output_file <- "Raw_Rain_Data.csv"
}



# --------------------------------------------------------------------------
#             Restart From Here
# --------------------------------------------------------------------------
# load the files
message("loading input data files")
weather_df_Final_comp <- read.csv(file=weather_dataset_in,header=TRUE,sep=",")
pollen_df_Final_comp <- read.csv(file=pollen_dataset_in,header=TRUE,sep=",")
rain_df_Final_comp <-read.csv(file=rain_dataset_in,header=TRUE,sep=",")
station_data_clean <- read.csv(file=station_dataset_in,header=TRUE,sep=",")

# --------------------------------------------------------------------
#   Match the postcodes/coordinates into the data from the station df
# --------------------------------------------------------------------
message("matching postcode and coordinate data to stations")
# transform into factors the df station
station_data_clean$Station <- 
  as.factor(station_data_clean$Station)
weather_df_Final_comp$Station <- 
  as.factor(weather_df_Final_comp$Station)
pollen_df_Final_comp$Station <- 
  as.factor(pollen_df_Final_comp$Station)
rain_df_Final_comp$Station <- 
  as.factor(rain_df_Final_comp$Station)

weather_df_Final_comp_post <- 
  merge(weather_df_Final_comp,station_data_clean,by="Station")
# we need to postpone the mergin with the pollen for later
#pollen_df_Final_comp_post <- 
  #merge(pollen_df_Final_comp,station_data_clean,by="Station")
rain_df_Final_comp_post <- 
  merge(rain_df_Final_comp,station_data_clean,by="Station")

#View(weather_df_Final_comp_post)
#View(weather_df_Final_comp)

#View(rain_df_Final_comp_post)
#View(rain_df_Final_comp)


#some checks ##############################################################
message("data checks:")
length(unique(weather_df_Final_comp_post$Date))#1095
length(unique(pollen_df_Final_comp$Date))#547
length(unique(rain_df_Final_comp_post$Date))#1095

length(unique(weather_df_Final_comp_post$Station))#519
length(unique(pollen_df_Final_comp$Station))#14 
length(unique(rain_df_Final_comp_post$Station))#2549

nrow(weather_df_Final_comp_post)/
  length(unique(weather_df_Final_comp_post$Station))#1095
nrow(pollen_df_Final_comp)/
  length(unique(pollen_df_Final_comp$Station))#547
nrow(rain_df_Final_comp_post)/
  length(unique(rain_df_Final_comp_post$Station))#1095
# There are some duplicates. How many rows should I have?
1095*2549 #2791155

# find the duplicates ######################################################
message("finding and removing duplicate")
#View(rain_df_Final_comp_post[duplicated(rain_df_Final_comp_post[,c(1,3)]),])
# elimnate the duplicates
#Indexes of the duplicate rows that will be removed: 
duplicate_indexes <- 
  which(duplicated(rain_df_Final_comp_post[c('Station', 'Date')]),) 
#duplicate_indexes 

#rain_df_Final_comp_post2 will contain unique dataset without the duplicates. 
rain_df_Final_comp_post2 <- rain_df_Final_comp_post[!duplicated(rain_df_Final_comp_post[c('Station', 'Date')]),] 
#View(rain_df_Final_comp_post2)
# check 
nrow(rain_df_Final_comp_post2)/
  length(unique(rain_df_Final_comp_post2$Station))#1095



message("extending pollen data to cover missing periods of data collection")
# Extend pollen_df_Final_comp by adding the missing periods of data collection ---
# After Imputation by Chained Equation populate the full date range for pollen ---
# Extend the pollen DF to include all possible combo of station and date -----
# create a DF with all the stations and dates
# list with all date for 3 years
date_list <- unique(weather_df_Final_comp$Date)
Station_rep <- c()
for(i in unique(pollen_df_Final_comp$Station)){
  L <- rep(i, length(date_list))
  Station_rep <- append(Station_rep, L)
}
Date_rep <- rep(date_list,length(unique(pollen_df_Final_comp$Station)))
StationDate_df <- data.frame(Station = Station_rep,
                             Date = Date_rep)

# before merging we need the Date column into the right format 2017-03-01 
StationDate_df$Date <- as.Date(StationDate_df$Date, 
                                          format = "%Y-%m-%d")
pollen_df_Final_comp$Date <- as.Date(pollen_df_Final_comp$Date, 
                             format = "%Y-%m-%d")

# merge the 2 dataframe to have all Stations and columns
pollen_df_Final_comp2 <- merge(pollen_df_Final_comp, StationDate_df, 
                          by = c("Station","Date"), all = TRUE)

# some checks 
nrow(pollen_df_Final_comp2)/
  length(unique(pollen_df_Final_comp2$Station)) # 1095
# ---------------------------------------------------------------------

message("mapping postcode and coordinates to the pollen stations")
# now we can map the postcode and coordinates to the pollen df
pollen_df_Final_comp_post <- 
  merge(pollen_df_Final_comp2,station_data_clean,by="Station")
#View(pollen_df_Final_comp_post)
#View(pollen_df_Final_comp2)

message("writing final data files")
# final DFs
#View(weather_df_Final_comp_post)
#View(pollen_df_Final_comp_post)
#View(rain_df_Final_comp_post2)
write.csv(weather_df_Final_comp_post, file = weather_output_file, row.names=F)
write.csv(pollen_df_Final_comp_post, file = pollen_output_file, row.names=F)
write.csv(rain_df_Final_comp_post2, file = rain_output_file, row.names=F)

