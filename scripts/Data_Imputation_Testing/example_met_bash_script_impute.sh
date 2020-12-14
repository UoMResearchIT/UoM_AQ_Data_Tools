#!/bin/bash

FILE_IN_TEST='inputs/Anns_for_Doug_tests/temp_extras-rel_hum-pressure-dewpoint_June2017.csv'
FILE_IN="inputs/temp_extras-rel_hum-pressure-dewpoint_midlands.csv"
OUT_DIR='met_postprocessing'
OUTFILE_SUFFIX='_midlands'
STATIONS_FILENAME="../../station_data/station_data_clean.csv"

DATE_RANGE="2017-01-01_00 2017-12-31_23"
MIN_YEARS=0.04  # 1
MIN_TEMP=-20
REFERENCE_NUM_YEARS=0.0625  # 3.5
REFERENCE_NUM_STATIONS=5
EXCLUDE_SITE_LIST='117'
VERBOSE=2

# EXAMPLE 1
ARGUMENTS_1=" --file_in ${FILE_IN_TEST}  --stations_filename ${STATIONS_FILENAME}
  --out_dir ${OUT_DIR} --outfile_suffix ${OUTFILE_SUFFIX} --date_range ${DATE_RANGE}
  --exclude_sites ${EXCLUDE_SITE_LIST} --min_temp ${MIN_TEMP}  --min_years ${MIN_YEARS}
  --ref_num_years ${REFERENCE_NUM_YEARS} --ref_num_stations ${REFERENCE_NUM_STATIONS}
  --impute_values --print_stats --verbose ${VERBOSE}"
  
  # TEST
  
FILE_IN_TEST='inputs/Anns_for_Doug_tests/Met_extracted_temp_extras-rel_hum-pressure-dewpoint_June2017.csv'
OUT_DIR_TEST='met_postprocessing_test_with_dougs'
OUTFILE_SUFFIX_TEST='_June2017_imputed'
DATE_RANGE_TEST="2017-01-01_00 2017-06-30_23"



FILE_IN="inputs/data_met/temp_rh_press_dewtemp_2016-2019.csv"
VERBOSE=1
MIN_YEARS=3
REFERENCE_NUM_YEARS=3.5
DATE_RANGE="2016-01-01_00 2019-12-31_23"
OUTFILE_SUFFIX='2016-2019_tests'



DATA_LOSS_POS='middle'
DATA_LOST='0.25'
SITES="513 3 1543 1033 1039 1046 23 24089 24090 18974 1055 32 61986"


ARGUMENTS_TEST="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}
  --file_in ${FILE_IN}  --stations_filename ${STATIONS_FILENAME}
  --out_dir ${OUT_DIR_TEST} --outfile_suffix ${OUTFILE_SUFFIX} --date_range ${DATE_RANGE}
  --exclude_sites ${EXCLUDE_SITE_LIST} --min_temp ${MIN_TEMP}  --min_years ${MIN_YEARS}
  --min_years_ref ${REFERENCE_NUM_YEARS} --ref_num_stations ${REFERENCE_NUM_STATIONS}
  --print_stats --verbose ${VERBOSE} --sites ${SITES}"


python Met_Temperature_RH_Pressure_impute_testing.py ${ARGUMENTS_TEST}
