#!/bin/bash

## Fixed Arguments. Adjust the verbosity (0=none) as required.
STATIONS_FILENAME="../../station_data/station_data_clean.csv"
EXCLUDE_SITE_LIST='117'
MIN_TEMP=-20
VERBOSE=2
REFERENCE_NUM_STATIONS=5
PRINT_STATS="--print_stats"
#PRINT_STATS=""

FIXED_ARGS="--stations_filename ${STATIONS_FILENAME} --min_temp ${MIN_TEMP} 
			--exclude_sites ${EXCLUDE_SITE_LIST} --ref_num_stations ${REFERENCE_NUM_STATIONS}
			--verbose ${VERBOSE} ${PRINT_STATS}"


## Comment / uncomment these lines to select the processing options you wish to use.
SCENARIO='REDUCED'
#SCENARIO='FULL'

#IMPUTATION=( 'IMPUTED' )
IMPUTATION=( 'NOT_IMPUTED' )



## Setting of scenario flags. Add to / edit these for your own scenarios.
if [[ $SCENARIO == 'FULL' ]]; then
	FILE_IN="inputs/data_met/temp_rh_press_dewtemp_2016-2019.csv"
	OUT_DIR='met_postprocessing'
	OUTFILE_SUFFIX='_2016-2019' # added to the imputation arguments

	DATE_RANGE="2016-01-01_00 2019-12-31_23"
	MIN_YEARS=2
	REFERENCE_NUM_YEARS=3.5

	SCENARIO_ARGS=" --file_in ${FILE_IN_TEST}  --out_dir ${OUT_DIR} --date_range ${DATE_RANGE}
					 --min_years ${MIN_YEARS} --ref_num_years ${REFERENCE_NUM_YEARS}"  
elif [[ $SCENARIO == 'REDUCED' ]]; then
	FILE_IN="inputs/temp_extras-rel_hum-pressure-dewpoint_midlands.csv"
	OUT_DIR='met_postprocessing'
	OUTFILE_SUFFIX='_midlands' # added to the imputation arguments
	
	DATE_RANGE="2017-01-01_00 2017-06-30_23"
	MIN_YEARS=0.04  # 1
	REFERENCE_NUM_YEARS=0.0625  # 3.5

	SCENARIO_ARGS=" --file_in ${FILE_IN_TEST}  --out_dir ${OUT_DIR} --date_range ${DATE_RANGE}
					 --min_years ${MIN_YEARS} --ref_num_years ${REFERENCE_NUM_YEARS}"
else
	echo "no valid scenario selected, printing help information."
	SCENARIO_ARGS="--help"
fi


## Setting of imputation flag.
if [[ ${IMPUTATION[0]} == 'IMPUTED' ]]; then
	OUTFILE_SUFFIX_IMPUTE="_with_imputation"
	IMPUTE_ARGS="--impute_values --outfile_suffix ${OUTFILE_SUFFIX}${OUTFILE_SUFFIX_IMPUTE}"	
elif [[ ${IMPUTATION[0]} == 'NOT_IMPUTED' ]]; then
	IMPUTE_ARGS="--no_impute_values --outfile_suffix ${OUTFILE_SUFFIX}"
else
	echo "no choice on imputation made, printing help information"
	IMPUTE_ARGS="--help"
fi



python Met_Temperature_RH_Pressure_data_processing.py ${FIXED_ARGS} ${SCENARIO_ARGS} ${IMPUTE_ARGS}
