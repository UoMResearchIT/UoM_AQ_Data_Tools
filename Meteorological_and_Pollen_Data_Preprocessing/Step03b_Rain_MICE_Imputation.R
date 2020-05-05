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
  package.args <- c("mice")
  loadPackages(package.args)
}

initialize()

message("finished initialization")


# set working dir
setwd("/mnt/iusers01/support/mbessdl2/scratch/Turing_AQ_Data_Prep")

# input file names
rain_dataset_in <- "rain_df_Max_Mean_Intermediate.csv"


# output file names
rain_dataset_out <- "rain_df_Final_comp.csv"




# --------------------------------------------------------------------------
#             Restart From Here
# --------------------------------------------------------------------------
# load the files
message("loading input data files")
rain_df_Final <-read.csv(file=rain_dataset_in,header=TRUE,sep=",")


# Apply MICE to impute missing values to individual files ------------ 
message("starting MICE calculations")
# (only for stations failure)
# Notes on MICE parameters:
#   m=5 refers to the number of imputed datasets. Five is the default value.
#   meth='pmm' refers to the imputation method. 
#   In this case we are using predictive mean matching as imputation method. 
#   Other imputation methods can be used, 
#   type methods(mice) for a list of the available imputation methods.

message("MICE: rain")
# rain (the rain df is too large, need to split: 2,910,046 obs, into 5)
# nrow(rain_df_Final)/5 #582009.2
# 1:582009
data3 <- rain_df_Final[1:582009,-2] 
tempData3 <- mice(data3,m=5,maxit=5,meth='pmm',seed=500)
summary(tempData3)
# 582010:1164019
data4 <- rain_df_Final[582010:1164019,-2] 
tempData4 <- mice(data4,m=5,maxit=5,meth='pmm',seed=500)
summary(tempData4)
# 1164020:1746029 
data5 <- rain_df_Final[1164020:1746029,-2] 
tempData5 <- mice(data5,m=5,maxit=5,meth='pmm',seed=500)
summary(tempData5)
# 1746030:2328039,
data6 <- rain_df_Final[1746030:2328039,-2] 
tempData6 <- mice(data6,m=5,maxit=5,meth='pmm',seed=500)
summary(tempData6)
# 2328040:2910046
data7 <- rain_df_Final[2328040:2910046,-2] 
tempData7 <- mice(data7,m=5,maxit=5,meth='pmm',seed=500)
summary(tempData7)
# complete the df
rain_df_Final_comp1 <- complete(tempData3,1)
rain_df_Final_comp2 <- complete(tempData4,1)
rain_df_Final_comp3 <- complete(tempData5,1)
rain_df_Final_comp4 <- complete(tempData6,1)
rain_df_Final_comp5 <- complete(tempData7,1)
# add the date to the DF
rain_df_Final_comp1$Date <- rain_df_Final[1:582009,]$Date
rain_df_Final_comp2$Date <- rain_df_Final[582010:1164019,]$Date
rain_df_Final_comp3$Date <- rain_df_Final[1164020:1746029,]$Date
rain_df_Final_comp4$Date <- rain_df_Final[1746030:2328039,]$Date
rain_df_Final_comp5$Date <- rain_df_Final[2328040:2910046,]$Date
# rbind the 5 dfs for rain
rain_df_Final_comp <- rbind(rain_df_Final_comp1,
                            rain_df_Final_comp2,
                            rain_df_Final_comp3,
                            rain_df_Final_comp4,
                            rain_df_Final_comp5)

message("writing final datafiles")
# temporarily save DF files #############################################
write.csv(rain_df_Final_comp, file = rain_dataset_out, row.names=F)
# --------------------------------------------------------------------


