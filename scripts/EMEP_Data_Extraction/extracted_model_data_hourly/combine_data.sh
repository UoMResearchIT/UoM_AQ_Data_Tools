#!/bin/bash

IN_START="emep_site_data_"
IN_TAIL=".csv"


DATES=( 'MarApr2016' 'MayJune2016' 'JulyAug2016' 'SeptOct2016' 'NovDec2016' \
		'JanFeb2017' 'MarApr2017' 'MayJune2017' 'JulyAug2017' 'SeptOct2017' 'NovDec2017' \
		'JanFeb2018' 'MarApr2018' 'MayJune2018' 'JulyAug2018' 'SeptOct2018' 'NovDec2018' \
		'JanFeb2019' 'MarApr2019' 'MayJune2019' 'JulyAug2019' 'SeptOct2019' 'NovDec2019' )

DATEVOIDS=( '2016-05-01' '2016-07-01' '2016-09-01' '2016-11-01' 2017-01-01' \
		'2017-03-01' '2017-05-01' '2017-07-01' '2017-09-01' '2017-11-01' 2018-01-01' \ 
		'2018-03-01' '2018-05-01' '2018-07-01' '2018-09-01' '2018-11-01' 2019-01-01' \
		'2019-03-01' '2019-05-01' '2019-07-01' '2019-09-01' '2019-11-01' 2020-01-01' )

OUT_FILE="emep_site_data_hourly_2016-2019.csv"
touch ${OUT_FILE}
echo -n > ${OUT_FILE}


first_file='emep_site_data_JanFeb2016.csv'

last_day=$( tail -1 $first_file | sed -n "s/^\([0-9]*-[0-9]*-[0-9]*\).*$/\1/p" )


grep -v ${last_day} $first_file >> ${OUT_FILE}


for DATE in ${DATES[@]}; do

	FILE_IN=${IN_START}${DATE}${IN_TAIL}

	last_day=$( tail -1 $FILE_IN | sed -n "s/^\([0-9]*-[0-9]*-[0-9]*\).*$/\1/p" )

	grep -v ${last_day} ${FILE_IN} | grep -v Date >> ${OUT_FILE} 

done
