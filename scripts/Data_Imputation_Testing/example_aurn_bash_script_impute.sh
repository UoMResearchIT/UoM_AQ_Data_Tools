#!/bin/bash

EMEP_FILE="emep_site_data_hourly_2016-2019.csv"
META_DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
META_DATA_FILENAME="AURN_metadata.RData"
SITES="ABD ABD7 ABD8 ARM6 AH"


# EXAMPLE 1
YEARS="2016 2017"
ARGUMENTS_1=" --metadata_url ${META_DATA_URL} --metadata_filename ${META_DATA_FILENAME} --years ${YEARS} --sites ${SITES} --verbose 2"

# EXAMPLE 2
YEARS="2017"
ARGUMENTS_2=" --metadata_url ${META_DATA_URL} --metadata_filename ${META_DATA_FILENAME} --years ${YEARS} --sites ${SITES} --verbose 2"




# For Dougs test

META_DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
META_DATA_FILENAME="AURN_metadata.RData"
OUT_DIR="../AURN_Data_Download/test_Ann"
#SITES="BARN BAR2 BAR3 BIL BIR BRAD BDMA CHS6"
SITES="LEED MID MY1 PT4"
SITES="LEED"
MINYEARS='1'
USEYEARS='3.8'
VERBOSE=3
DATE_RANGE="2016-01-01_00 2019-12-31_23"
INFILE_SUFFIX='2016-2017'
OUTFILE_SUFFIX='2016-2017_imputed'
DATA_LOSS_POS='middle'
DATA_LOST='0.25'

ARGUMENTS_1="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST} --metadata_url ${META_DATA_URL}
             --date_range ${DATE_RANGE} --min_years ${MINYEARS} --min_years_ref ${USEYEARS} 
             --metadata_filename ${META_DATA_FILENAME} --sites ${SITES} -o ${OUT_DIR} 
             -s ${INFILE_SUFFIX} --verbose ${VERBOSE}"

#echo ${ARGUMENTS_1}

python AURN_impute_testing.py ${ARGUMENTS_1}
