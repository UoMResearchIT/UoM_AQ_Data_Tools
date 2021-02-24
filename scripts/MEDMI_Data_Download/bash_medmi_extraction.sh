#!/bin/bash


#SCENARIO="FULL"
SCENARIO="REDUCED"

if [[ $SCENARIO == 'FULL' ]]; then

	LATITUDES='48 60'
	LONGITUDES='-11 3'
	VERBOSE=0
	MEASUREMENTS='temperature rain wind pollen'
	DATE_RANGE='2016-01-01_0 2019-12-31_23'
	OUTDIR="full_data"
	OUTFILE_SUFFIX="full"
	EXTRA_DATA_FLAG="--extra_measurements"

elif [[ $SCENARIO == 'REDUCED' ]]; then

	LATITUDES='55 56'
	LONGITUDES='-5 -3'
	VERBOSE=1
	MEASUREMENTS='temperature rain wind pollen'
	DATE_RANGE='2017-06-01_0 2017-06-30_23'
	OUTDIR="reduced_data"
	OUTFILE_SUFFIX="reduced"
	EXTRA_DATA_FLAG="--extra_measurements"

fi


ARGUMENTS_DATASET="--measurements ${MEASUREMENTS} ${EXTRA_DATA_FLAG}"
ARGUMENTS_GEOTEMP="--date_range ${DATE_RANGE} --latitude_range ${LATITUDES} --longitude_range ${LONGITUDES}"
ARGUMENTS_FILE="--outdir ${OUTDIR} --outfile_suffix ${OUTFILE_SUFFIX} -v ${VERBOSE}"


python met_extraction_script.py ${ARGUMENTS_DATASET} ${ARGUMENTS_GEOTEMP} ${ARGUMENTS_FILE}
