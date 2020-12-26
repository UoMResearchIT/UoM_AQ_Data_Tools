#!/bin/bash 


SENSOR="./AURN_metadata.RData"
SENSORTYPE="RDATA"

WRF="../../../WRF_UK/2017/wrfout_d01_2017-05-25_00:00:00"

EMEP="../../../EMEP_UK/Base_hourInst_JanFeb2016.nc"


# Extracting Daily Mean / Max Data
OUT="./extracted_model_data_daily/emep_site_data_JanFeb2016.csv"
DAILY_FLAG="daily_means"
ARGUMENTS=" --sensor_file ${SENSOR} --sensor_file_type ${SENSORTYPE} --wrf_file ${WRF} --emep_file ${EMEP} --out_file ${OUT} --${DAILY_FLAG}"
python emep_extract.py ${ARGUMENTS}



# Extracting the original, hourly, data
OUT="./extracted_model_data_hourly/emep_site_data_JanFeb2016.csv"
DAILY_FLAG="no_daily_means"
ARGUMENTS=" --sensor_file ${SENSOR} --sensor_file_type ${SENSORTYPE} --wrf_file ${WRF} --emep_file ${EMEP} --out_file ${OUT} --${DAILY_FLAG}"
python emep_extract.py ${ARGUMENTS}
