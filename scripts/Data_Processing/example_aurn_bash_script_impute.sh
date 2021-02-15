#!/bin/bash

EMEP_FILE="emep_site_data_hourly_2016-2019.csv"
META_DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
META_DATA_FILENAME="AURN_metadata.RData"
VERBOSE=2

FIXED_ARGS="--metadata_url ${META_DATA_URL} --metadata_filename ${META_DATA_FILENAME} --verbose ${VERBOSE}"


SCENARIO='REDUCED'
#SCENARIO='FULL'

#IMPUTATION=( 'IMPUTED' 'WITH_EMEP' )
#IMPUTATION=( 'IMPUTED' 'NO_EMEP' )
IMPUTATION=( 'NOT_IMPUTED' )



if [[ $SCENARIO == 'FULL' ]]; then

	OUT_DIR="../AURN_Data_Download/AURN_datasets"
	INFILE_SUFFIX='full_2016-2019'
	OUTFILE_SUFFIX='full_2016-2019'
	DATE_RANGE="2016-01-01_00 2019-12-31_23"
	MINYEARS='2'
	REFYEARS='3.5'

	SCENARIO_ARGS=" --date_range ${DATE_RANGE} --outdir_name ${OUT_DIR} --outfile_suffix ${OUTFILE_SUFFIX}
	                --infile_suffix ${INFILE_SUFFIX} --min_years ${MINYEARS} --min_years_ref ${REFYEARS} "

elif [[ $SCENARIO == 'REDUCED' ]]; then

	OUT_DIR="../AURN_Data_Download/AURN_datasets"
	SITES="BARN BAR2 BAR3 BIL BIR BRAD BDMA CHS6 CHLG CHS7 DYAG DCST FEA HSAW"
	INFILE_SUFFIX='reduced_2016-2017'
	OUTFILE_SUFFIX='reduced_2016-2017'
	DATE_RANGE="2016-01-01_00 2017-12-31_23"
	MINYEARS='1'
	REFYEARS='1'

	SCENARIO_ARGS=" --date_range ${DATE_RANGE} --outdir_name ${OUT_DIR} --outfile_suffix ${OUTFILE_SUFFIX}
	                --infile_suffix ${INFILE_SUFFIX} --min_years ${MINYEARS} --min_years_ref ${REFYEARS} 
	                --sites ${SITES} "

else
	echo "no valid scenario selected, printing help information."
	SCENARIO_ARGS="--help"
fi


if [[ ${IMPUTATION[0]} == 'IMPUTED' ]]; then

	if [[ ${IMPUTATION[1]} == 'WITH_EMEP' ]]; then
		
		EMEP_FILENAME="../EMEP_Data_Extraction/extracted_model_data_hourly/emep_site_data_hourly_2016-2019.csv"
		
		IMPUTE_ARGS="--impute_values --emep_filename ${EMEP_FILENAME}"
		
	else 
		IMPUTE_ARGS="--impute_values"	
	fi

elif [[ ${IMPUTATION[0]} == 'NOT_IMPUTED' ]]; then

	IMPUTE_ARGS="--no_impute_values"

else
	echo "no choice on imputation made, printing help information"
	IMPUTE_ARGS="--help"
fi



python AURN_processing.py ${FIXED_ARGS} ${SCENARIO_ARGS} ${IMPUTE_ARGS}
