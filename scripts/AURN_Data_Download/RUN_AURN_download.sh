#!/bin/bash

META_DATA_URL="https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
META_DATA_FILENAME="AURN_metadata.RData"
VERBOSE=2


FIXED_ARGS="--metadata_url ${META_DATA_URL} --metadata_filename ${META_DATA_FILENAME} --verbose ${VERBOSE}"

SCENARIO='REDUCED'
#SCENARIO='FULL'

if [[ $SCENARIO == 'FULL' ]]; then

	OUT_DIR="AURN_datasets"
	YEARS="2016 2017 2018 2019"
	OUTFILE_SUFFIX='full_2016-2019'

	SCENARIO_ARGS=" --years ${YEARS} -o ${OUT_DIR} -s ${OUTFILE_SUFFIX}"


elif [[ $SCENARIO == 'REDUCED' ]]; then

	OUT_DIR="AURN_datasets"
	YEARS="2016 2017"
	SITES="BARN BAR2 BAR3 BIL BIR BRAD BDMA CHS6 CHLG CHS7 DYAG DCST FEA HSAW"
	OUTFILE_SUFFIX='reduced_2016-2017'

	SCENARIO_ARGS="--years ${YEARS} -o ${OUT_DIR} -s ${OUTFILE_SUFFIX} --sites ${SITES}"

else
	echo "no valid scenario selected, printing help information."
	SCENARIO_ARGS="--help"
fi


echo 'flags used:' ${FIXED_ARGS} ${SCENARIO_ARGS}
python AURN_download.py ${FIXED_ARGS} ${SCENARIO_ARGS}
