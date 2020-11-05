#!/bin/bash

EMEP_FILE="emep_site_data_hourly_2016-2019.csv"
META_DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
META_DATA_FILENAME="AURN_metadata.RData"
SITES="ABD ABD7 ABD8 ARM6 AH"


# EXAMPLE 1
YEARS="2016 2017"
ARGUMENTS_1=" --metadata_url ${META_DATA_URL} --metadata_filename ${META_DATA_FILENAME} --years ${YEARS} --sites ${SITES}"

# EXAMPLE 2
YEARS="2017"
ARGUMENTS_2=" --metadata_url ${META_DATA_URL} --metadata_filename ${META_DATA_FILENAME} --years ${YEARS} --sites ${SITES} --verbose 2"


python AURN_download.py ${ARGUMENTS_2}


