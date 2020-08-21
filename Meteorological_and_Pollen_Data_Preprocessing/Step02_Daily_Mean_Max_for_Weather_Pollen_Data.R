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

initialize()

message("finished initialization")


pollen_file<-"midas_pollen_drnl_ob.clean.csv"
weather_17_file<-"midas_weather_hrly_ob__2017.clean.csv"
weather_18_file<-"midas_weather_hrly_ob__2018.clean.csv"
weather_19_file<-"midas_weather_hrly_ob__2019.clean.csv"
rain_17_file<-"midas_rain_drnl_ob__1__2017.clean.csv"
rain_18_file<-"midas_rain_drnl_ob__1__2018.clean.csv"
rain_19_file<-"midas_rain_drnl_ob__1__2019.clean.csv"



weather_out_file<-"weather_df_Max_Mean_Intermediate.csv"
pollen_out_file<-"pollen_df_Max_Mean_Intermediate.csv"
rain_out_file<-"rain_df_Max_Mean_Intermediate.csv"


# ----------------------------------------------------------------------
#                 Load and clean the datasets  
# ----------------------------------------------------------------------
# set working dir
setwd("/Users/mbessdl2/work/manchester/Turing_Project_AQ/Britain_Breathing_Data/MIDAS_Datasets_Processing/Manuele_Reani_Raw_Data/test_processes")
# load tha data 
message("loading pollen data")
pollen_data <- read.csv(file=pollen_file,header=TRUE,sep=",")

message("loading and merging weather data")
weather_data17 <- read.csv(file=weather_17_file,header=TRUE,sep=",")
weather_data18 <- read.csv(file=weather_18_file,header=TRUE,sep=",")
weather_data19 <- read.csv(file=weather_19_file,header=TRUE,sep=",")
names(weather_data18)<-names(weather_data17)
names(weather_data19)<-names(weather_data17)
weather_data <- rbind(weather_data17, weather_data18, weather_data19)

message("loading and merging rain data")
rain_data17 <- read.csv(file=rain_17_file,header=TRUE,sep=",")
rain_data18 <- read.csv(file=rain_18_file,header=TRUE,sep=",")
rain_data19 <- read.csv(file=rain_19_file,header=TRUE,sep=",")
names(rain_data18)<-names(rain_data17)
names(rain_data19)<-names(rain_data17)
rain_data <- rbind(rain_data17, rain_data18, rain_data19) 

message("dropping non-required weather data")
# only keep relevant features for weather
weather_data <- weather_data[,c(1,2,29,32,33)]
message("labelling weather and rain data")
# name these features 
names(weather_data) <- 
  c("Station", "Date.precise", "Temperature", "Humidity", "Pressure")
# name the features of Rain df 
names(rain_data) <- 
  c("Station", "Date.precise", "Daily.rain") 

message("deleting empty records")
# eliminate empty record
#pollen_data <- pollen_data[-c(10000, 20000, 22300, 22301),]
#weather_data <- weather_data[-2950071, -6026440, -9748861]
#rain_data <- rain_data[-2564857, -1735098, -871592] 

message("removing timestamp from date information")
# add the Date column and remove the Date.precise 
# weather 2017-01-01 00:00:00
weather_data$Date <- as.Date(weather_data$Date.precise, 
                             format = "%Y-%m-%d")
weather_data <- weather_data[,-2]
# rain 2017-01-01 00:00:00
rain_data$Date <- as.Date(rain_data$Date.precise, 
                             format = "%Y-%m-%d")
rain_data <- rain_data[,-2]

# pollen 20/01/2010 00:00
pollen_data$Date <- as.Date(pollen_data$Date.precise, 
                             format = "%d/%m/%Y")
pollen_data <- pollen_data[,-2]

message("selecting data for 2017-2019 only")
# only keep data from 2017 - 2019 for the 3 datasets 
weather_data <- weather_data[weather_data$Date >= "2017-01-01" & 
               weather_data$Date <= "2019-12-31",]
rain_data <- rain_data[rain_data$Date >= "2017-01-01" & 
            rain_data$Date <= "2019-12-31",]
pollen_data <- pollen_data[pollen_data$Date >= "2017-01-01" & 
              pollen_data$Date <= "2019-12-31",]

message("deleting excess empty rows")
# eliminate the rows with all NAs at the bottom of the df
weather_data <- weather_data[!(is.na(weather_data$Station)),]
pollen_data <- pollen_data[!(is.na(pollen_data$Station)),]
rain_data <- rain_data[!(is.na(rain_data$Station)),]

# ---------------------------------------------------------------------
#           WEATHER DATA aggregate mean/max and dataset extension 
# ---------------------------------------------------------------------

message("calculating daily mean and max values")
# average and max for each day on weather -----------------------------
mean_weather <- aggregate(weather_data[,2:4],
                          by=list(weather_data$Station, weather_data$Date),
                          FUN=mean, na.rm = TRUE)
max_weather <- aggregate(weather_data[,2:4],
                          by=list(weather_data$Station, weather_data$Date),
                          FUN=max, na.rm = TRUE)
names(mean_weather) <- c("Station", "Date", 
                         "Temperature.mean", "Humidity.mean", "Pressure.mean")
names(max_weather) <-  c("Station", "Date", 
                         "Temperature.max", "Humidity.max", "Pressure.max")

message("replacing NaN and Infinity with NA")
# replace NaN with NA
is.na(mean_weather) <- sapply(mean_weather, is.nan)
# replace infinity with NA
is.na(max_weather) <- sapply(max_weather, is.infinite)

message("Rounding to 1st decimal place")
# round the mean df to the 1st decimal 
mean_weather <- mean_weather %>% mutate_at(vars(Temperature.mean, 
                                   Humidity.mean,
                                   Pressure.mean), 
                              funs(round(., 1)))

message("Merging Weather Mean and Max Datasets")
# merge the two data frames by Station and Date
weather_df <- merge(mean_weather,max_weather,by=c("Station", "Date"))
# ---------------------------------------------------------------------

message("Extending weather timeseries for all stations")
# Extend the DF to include all possible combo of station and date -----
# create a DF with all the stations and dates 
Station_rep <- c()
for(i in unique(weather_df$Station)){
  L <- rep(i, length(unique(weather_df$Date)))
  Station_rep <- append(Station_rep, L)
}
Date_rep <- rep(unique(weather_df$Date),length(unique(weather_df$Station)))
StationDate_df <- data.frame(Station = Station_rep,
                             Date = Date_rep)
# merge the 2 dataframe to have all Stations and columns
weather_df_Final <- merge(weather_df, StationDate_df, 
                 by = c("Station","Date"), all = TRUE)
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
#                 POLLEN DATA dataset extension 
# ---------------------------------------------------------------------

message("Extending pollen timeseries for all stations")
# Extend the DF to include all possible combo of station and date -----
# create a DF with all the stations and dates
Station_rep <- c()
for(i in unique(pollen_data$Station)){
  L <- rep(i, length(unique(pollen_data$Date)))
  Station_rep <- append(Station_rep, L)
}
Date_rep <- rep(unique(pollen_data$Date),length(unique(pollen_data$Station)))
StationDate_df <- data.frame(Station = Station_rep,
                             Date = Date_rep)
# merge the 2 dataframe to have all Stations and columns
pollen_df_Final <- merge(pollen_data, StationDate_df, 
                         by = c("Station","Date"), all = TRUE)
# ------------------------------------------------------------------------

# ---------------------------------------------------------------------
#                 RAIN DATA dataset extension 
# ---------------------------------------------------------------------

message("Extending rain timeseries for all stations")
# Extend the DF to include all possible combo of station and date -----
# create a DF with all the stations and dates
Station_rep <- c()
for(i in unique(rain_data$Station)){
  L <- rep(i, length(unique(rain_data$Date)))
  Station_rep <- append(Station_rep, L)
}
Date_rep <- rep(unique(rain_data$Date),length(unique(rain_data$Station)))
StationDate_df <- data.frame(Station = Station_rep,
                             Date = Date_rep)
# merge the 2 dataframe to have all Stations and columns
rain_df_Final <- merge(rain_data, StationDate_df, 
                         by = c("Station","Date"), all = TRUE)
# ----------------------------------------------------------------------
 
message("Saving Intermediate Data Files")
# temporarily save DF files #############################################
write.csv(weather_df_Final, file = weather_out_file, row.names=F)
write.csv(pollen_df_Final, file = pollen_out_file, row.names=F)
write.csv(rain_df_Final, file = rain_out_file, row.names=F)
# --------------------------------------------------------------------





