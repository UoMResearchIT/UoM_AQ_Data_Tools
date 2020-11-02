#!/bin/bash

# FILE_IN = 'inputs/temp_extras-rel_hum-pressure-dewpoint_June2017.csv'
FILE_IN="inputs/temp_extras-rel_hum-pressure-dewpoint_midlands.csv"
OUT_DIR='met_postprocessing'
OUTFILE_SUFFIX='_midlands'
STATIONS_FILENAME="../station_data/station_data_clean.csv"

DATE_RANGE='2017-01-01_00 2017-12-31_23'
MIN_YEARS=0.04  # 1
MIN_TEMP=-20
REFERENCE_NUM_YEARS=0.0625  # 3.5
REFERENCE_NUM_STATIONS=5
EXCLUDE_SITE_LIST=117
VERBOSE=2

# EXAMPLE 1
ARGUMENTS_1=" --file_in ${FILE_IN}  --stations_filename ${STATIONS_FILENAME}
  --out_dir ${OUT_DIR} --outfile_suffix ${OUTFILE_SUFFIX} --date_range ${DATE_RANGE}
  --exclude_sites ${EXCLUDE_SITE_LIST} --min_temp ${MIN_TEMP}  --min_years ${MIN_YEARS}
  --ref_num_years ${REFERENCE_NUM_YEARS} --ref_num_stations ${REFERENCE_NUM_STATIONS}
  --impute_values --print_stats --verbose ${VERBOSE}"

python Met_Temperature_RH_Pressure_data_processing.py ${ARGUMENTS_1}


