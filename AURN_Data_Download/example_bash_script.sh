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
OUT_DIR="test_Ann"
SITES="BARN BAR2 BAR3 BIL BIR BRAD BDMA CHS6 CHLG CHS7 DYAG DCST FEA HSAW HM HULL HUL2 HULR IMGM LB LEED LED6 LDS LIN3 LINC
MID NEWC NCA3 REDC ROTH SCUN SCN2 SHBR SHE2 SHDG SHE SOTR EAGL YARM SUND SUN2 SUNR YK10 YK11"
MINYEARS='1'
USEYEARS='1'
YEARS="2016 2017"
VERBOSE=3
OUTFILE_SUFFIX='2016-2017_no-impute' 

ARGUMENTS_1="--no_impute_values --metadata_url ${META_DATA_URL} --min_years ${MINYEARS} --min_years_ref ${USEYEARS} --metadata_filename ${META_DATA_FILENAME} --years ${YEARS} --sites ${SITES} -o ${OUT_DIR} -s ${OUTFILE_SUFFIX}  --verbose ${VERBOSE}"

python AURN_download.py ${ARGUMENTS_1}
