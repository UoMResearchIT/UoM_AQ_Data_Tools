#!/bin/bash

EMEP_FILE="emep_site_data_hourly_2016-2019.csv"
META_DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
META_DATA_FILENAME="AURN_metadata.RData"
DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/"
YEARS="2016 2017"
SITES="ABD ABD7 ABD8 ARM6 AH"


# EXAMPLE 1

ARGUMENTS_1=" --meta_data_url ${META_DATA_URL} --meta_data_filename ${META_DATA_FILENAME} --data_url ${DATA_URL} --years ${YEARS} --sites ${SITES}"

# EXAMPLE 2
#YEARS="2017"
ARGUMENTS_2=" --meta_data_url ${META_DATA_URL} --meta_data_filename ${META_DATA_FILENAME} --data_url ${DATA_URL} --years ${YEARS}"


python AURN_download.py ${ARGUMENTS_2}


