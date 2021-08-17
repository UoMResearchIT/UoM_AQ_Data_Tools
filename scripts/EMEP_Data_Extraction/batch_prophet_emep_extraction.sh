#!/bin/bash --login
#$ -N emep_extraction            
#$ -cwd
#$ -t 1-1


source ~/bin/conda_activate.sh 
conda activate wrf_emep_work


INDEX=$((SGE_TASK_ID-1))


DATES=( 'JanFeb2016' 'MarApr2016' 'MayJune2016' 'JulyAug2016' 'SeptOct2016' 'NovDec2016' \
		'JanFeb2017' 'MarApr2017' 'MayJune2017' 'JulyAug2017' 'SeptOct2017' 'NovDec2017' \
		'JanFeb2018' 'MarApr2018' 'MayJune2018' 'JulyAug2018' 'SeptOct2018' 'NovDec2018' \
		'JanFeb2019' 'MarApr2019' 'MayJune2019' 'JulyAug2019' 'SeptOct2019' 'NovDec2019' )

DATES=( 'MarAprMay2020' )

SENSOR="./sensors_detailed.csv"

WRF="../WRF_UK/2017/wrfout_d01_2017-05-25_00:00:00"

EMEP_PATH="../EMEP_UK/"
EMEP_START="Base_hourInst_"
EMEP_TAIL=".nc"

OUT_PATH="./extracted_model_data/"
OUT_START="emep_prophet_NO2_data_"
OUT_TAIL=".csv"


# thread specific data
DATE=${DATES[$INDEX]}
EMEP=${EMEP_PATH}${EMEP_START}${DATE}${EMEP_TAIL}
OUT=${OUT_PATH}${OUT_START}${DATE}${OUT_TAIL}

ARGUMENTS=" --sensor_file ${SENSOR} --wrf_file ${WRF} --emep_file ${EMEP} --out_file ${OUT}"

echo "arguments are: ${ARGUMENTS}"
python emep_extract_prophet.py ${ARGUMENTS}

echo "finished extraction task"
