#!/bin/bash --login
#$ -N wrf_extraction            
#$ -cwd
#$ -t 1-12


source ~/bin/conda_activate.sh 
conda activate wrf_emep_work


INDEX=$((SGE_TASK_ID-1))


MONTHS=( '01' '02' '03' '04' '05' '06' \
		'07' '08' '09' '10' '11' '12' )



SENSOR="./sensors_detailed.csv"

WRF_PATH_BASE="../WRF_UK/"


OUT_PATH="./test_data/"


pwd

# thread specific data
MONTH=${MONTHS[$INDEX]}

YEARS=( '2016' '2017' '2018' '2019' )

for YEAR in ${YEARS[@]}; do

	DATE=${YEAR}-${MONTH}
	WRF_PATH=${WRF_PATH_BASE}${YEAR}


	ARGUMENTS=" --sensor_file ${SENSOR} --wrf_file_path ${WRF_PATH} --date_range ${DATE} --out_file_path ${OUT_PATH}"

	echo "arguments are: ${ARGUMENTS}"
	python wrf_extract_pblh.py ${ARGUMENTS}


done

echo "finished extraction task"
