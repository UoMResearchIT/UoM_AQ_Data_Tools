#!/bin/bash

OUTFILE_PREFIX = 'dataxxx_'
DATE_RANGE = ['2016-1-1 0','2019-12-31 23']
DATES_STRING = '2016-2019'
LATITUDES = [48,60]
LONGITUDES = [-11,3]
VERBOSE = 0


# EXAMPLE 1
ARGUMENTS_1=" --outfile_prefix ${OUTFILE_PREFIX} --date_range ${DATE_RANGE} --dates_string ${DATE_STRING}"


python met_extraction_script.py ${ARGUMENTS_1}


