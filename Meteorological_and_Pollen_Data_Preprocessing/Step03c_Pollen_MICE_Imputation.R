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
pollen_dataset_in <- "pollen_df_Max_Mean_Intermediate.csv"


# output file names
pollen_dataset_out <- "pollen_df_Final_comp.csv"




# --------------------------------------------------------------------------
#             Restart From Here
# --------------------------------------------------------------------------
# load the files
message("loading input data files")
pollen_df_Final <- read.csv(file=pollen_dataset_in,header=TRUE,sep=",")


# Apply MICE to impute missing values to individual files ------------ 
message("starting MICE calculations")
# (only for stations failure)
# Notes on MICE parameters:
#   m=5 refers to the number of imputed datasets. Five is the default value.
#   meth='pmm' refers to the imputation method. 
#   In this case we are using predictive mean matching as imputation method. 
#   Other imputation methods can be used, 
#   type methods(mice) for a list of the available imputation methods.

message("MICE: pollen")
# pollen
data2 <- pollen_df_Final[,-2]
tempData2 <- mice(data2,m=5,maxit=5,meth='pmm',seed=500)
summary(tempData2)
pollen_df_Final_comp <- complete(tempData2,1)
pollen_df_Final_comp$Date <- pollen_df_Final$Date


message("writing final datafiles")
# temporarily save DF files #############################################
write.csv(pollen_df_Final_comp, file = pollen_dataset_out, row.names=F)
# --------------------------------------------------------------------


