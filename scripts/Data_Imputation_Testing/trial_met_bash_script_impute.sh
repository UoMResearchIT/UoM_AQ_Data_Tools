#!/bin/bash

STATIONS_FILENAME="../../station_data/station_data_clean.csv"
FILE_IN="inputs/data_met/temp_rh_press_dewtemp_2016-2019.csv"
DATE_RANGE="2016-01-01_00 2019-12-31_23"
OUTFILE_SUFFIX='2016-2019_tests'
OUT_DIR='met_postprocessing_test_with_dougs'
MIN_YEARS=3
REFERENCE_NUM_YEARS=3.5
REFERENCE_NUM_STATIONS=5
EXCLUDE_SITE_LIST='117'
MIN_TEMP=-20

VERBOSE=1

ARGUMENTS_FIXED="
  --file_in ${FILE_IN}  --stations_filename ${STATIONS_FILENAME}
  --out_dir ${OUT_DIR} --outfile_suffix ${OUTFILE_SUFFIX} --date_range ${DATE_RANGE}
  --exclude_sites ${EXCLUDE_SITE_LIST} --min_temp ${MIN_TEMP}  --min_years ${MIN_YEARS}
  --min_years_ref ${REFERENCE_NUM_YEARS} --ref_num_stations ${REFERENCE_NUM_STATIONS}
  --print_stats --verbose ${VERBOSE} 
"


run_analysis () 
{
DO_NOT_USE_ALL_SITES="513 3 1543 1033 1039 1046 23 24089 24090 18974 1055 32 1060 1575 44 556
 1070 1074 56370 52 1076 54 1078 1083 24125 67 583 56905 56906 56907 56908 56909
 79 595 62041 605 613 105 113 1137 1144 1145 643 132 56963 137 1161 1171 150 1180
 669 161 1190 56486 1198 177 692 1209 1215 17089 17090 17091 708 709 17094 17097
 1226 17098 17099 16589 17101 17102 24275 212 16596 726 16611 743 744 1255 57063
 235 19187 19188 19206 775 779 268 1302 17176 795 1319 811 1336 315 19260 30523
 1346 1352 842 847 16725 346 862 869 1383 1386 876 57199 1393 370 1395 373 886
 888 889 23417 381 384 386 1415 393 395 908 405 409 1435 30620 17309 17314
 24996 421 1448 1450 440 1467 17344 461 471 18903 987 30690 1007 498 1534 1023"

SITES_A="19188 726 744 842 847"
STAT_DIR_SITES_A="met_stats/test_stats_sitesA"
ARGUMENTS_SITES_A="--sites ${SITES_A} --statdir_name ${STAT_DIR_SITES_A}${SCEN_TAIL}"

SITES_B="52 743 17309 1534 847"
STAT_DIR_SITES_B="met_stats/test_stats_sitesB"
ARGUMENTS_SITES_B="--sites ${SITES_B} --statdir_name ${STAT_DIR_SITES_B}${SCEN_TAIL}"

SITES_C="161 1074 23 708 1302"
STAT_DIR_SITES_C="met_stats/test_stats_sitesC"
ARGUMENTS_SITES_C="--sites ${SITES_C} --statdir_name ${STAT_DIR_SITES_C}${SCEN_TAIL}"


python Met_Temperature_RH_Pressure_impute_testing.py ${ARGUMENTS_FIXED} ${ARGUMENTS_DATA} ${ARGUMENTS_SITES_A}
python Met_Temperature_RH_Pressure_impute_testing.py ${ARGUMENTS_FIXED} ${ARGUMENTS_DATA} ${ARGUMENTS_SITES_B}
python Met_Temperature_RH_Pressure_impute_testing.py ${ARGUMENTS_FIXED} ${ARGUMENTS_DATA} ${ARGUMENTS_SITES_C}
}



### base setup
DATA_LOSS_POS='random'
DATA_LOST='0.5'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_random_50"
#run_analysis

DATA_LOSS_POS='random'
DATA_LOST='0.75'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_random_75"
#run_analysis

DATA_LOSS_POS='random'
DATA_LOST='0.25'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_random_25"
#run_analysis

DATA_LOSS_POS='start'
DATA_LOST='0.5'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_startloss_50"
#run_analysis

DATA_LOSS_POS='start'
DATA_LOST='0.25'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_startloss_25"
#run_analysis

DATA_LOSS_POS='start'
DATA_LOST='0.75'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_startloss_75"
#run_analysis

DATA_LOSS_POS='end'
DATA_LOST='0.5'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_endloss_50"
run_analysis

DATA_LOSS_POS='end'
DATA_LOST='0.25'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_endloss_25"
#run_analysis

DATA_LOSS_POS='end'
DATA_LOST='0.75'
ARGUMENTS_DATA="--data_loss_position ${DATA_LOSS_POS} --data_lost ${DATA_LOST}"
SCEN_TAIL="_endloss_75"
#run_analysis

