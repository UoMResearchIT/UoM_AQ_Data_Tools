#!/bin/bash

## Fixed arguments. Adjust the verbosity (0=none) as required.
META_DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
META_DATA_FILENAME="AURN_metadata.RData"
VERBOSE=2

FIXED_ARGS="--metadata_url ${META_DATA_URL} --metadata_filename ${META_DATA_FILENAME} --verbose ${VERBOSE}"


## Comment / uncomment these lines to select the processing options you wish to use.
SCENARIO='REDUCED'
#SCENARIO='FULL'

#IMPUTATION=( 'IMPUTED' 'WITH_EMEP' )
#IMPUTATION=( 'IMPUTED' 'NO_EMEP' )
IMPUTATION=( 'NOT_IMPUTED' )


## Setting of scenario flags. Add to / edit these for your own scenarios.
if [[ $SCENARIO == 'FULL' ]]; then
	OUT_DIR="../AURN_Data_Download/AURN_datasets"
	INFILE_SUFFIX='full_2016-2019'
	OUTFILE_SUFFIX='full_2016-2019' # added to the imputation arguments
	DATE_RANGE="2016-01-01_00 2019-12-31_23"
	MINYEARS='2'
	REFYEARS='3.5'

	SCENARIO_ARGS=" --date_range ${DATE_RANGE} --outdir_name ${OUT_DIR} 
	                --infile_suffix ${INFILE_SUFFIX} --min_years ${MINYEARS} --min_years_ref ${REFYEARS} "
elif [[ $SCENARIO == 'REDUCED' ]]; then
	OUT_DIR="../AURN_Data_Download/AURN_datasets"
	SITES="BARN BAR2 BAR3 BIL BIR BRAD BDMA CHS6 CHLG CHS7 DYAG DCST FEA HSAW"
	INFILE_SUFFIX='reduced_2016-2017'
	OUTFILE_SUFFIX='reduced_2016-2017' # added to the imputation arguments
	DATE_RANGE="2016-01-01_00 2017-12-31_23"
	MINYEARS='1'
	REFYEARS='1'

	SCENARIO_ARGS=" --date_range ${DATE_RANGE} --outdir_name ${OUT_DIR}
	                --infile_suffix ${INFILE_SUFFIX} --min_years ${MINYEARS} --min_years_ref ${REFYEARS} 
	                --sites ${SITES} "
else
	echo "no valid scenario selected, printing help information."
	SCENARIO_ARGS="--help"
fi

## Setting of imputation flags, including EMEP data location and out file suffix
if [[ ${IMPUTATION[0]} == 'IMPUTED' ]]; then
	if [[ ${IMPUTATION[1]} == 'WITH_EMEP' ]]; then
		EMEP_FILENAME="../EMEP_Data_Extraction/extracted_model_data_hourly/emep_site_data_hourly_2016-2019.csv"
		OUTFILE_SUFFIX_IMPUTE="_with_imputation_inc_emep"
		
		IMPUTE_ARGS="--impute_values --emep_filename ${EMEP_FILENAME} 
						--outfile_suffix ${OUTFILE_SUFFIX}${OUTFILE_SUFFIX_IMPUTE}"
	else 
		OUTFILE_SUFFIX_IMPUTE="_with_imputation_inc_emep"

		IMPUTE_ARGS="--impute_values --outfile_suffix ${OUTFILE_SUFFIX}${OUTFILE_SUFFIX_IMPUTE}"	
	fi
elif [[ ${IMPUTATION[0]} == 'NOT_IMPUTED' ]]; then
	IMPUTE_ARGS="--no_impute_values --outfile_suffix ${OUTFILE_SUFFIX}"
else
	echo "no choice on imputation made, printing help information"
	IMPUTE_ARGS="--help"
fi


## Running of AURN processing script.
python AURN_processing.py ${FIXED_ARGS} ${SCENARIO_ARGS} ${IMPUTE_ARGS}
