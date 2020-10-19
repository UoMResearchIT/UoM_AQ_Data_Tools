#!/bin/bash

OUTDIR_PREFIX = 'dataxxx_'
OUTFILE_SUFFIX = '2017-01--2017-06'
DATE_RANGE = '2017-1-1 0 2017-06-30 23'
LATITUDES = 53 55
LONGITUDES = -5 -3
VERBOSE = 0


# EXAMPLE 1
ARGUMENTS_1=" --outdir_prefix ${OUTDIR_PREFIX} --outfile_suffix ${OUTFILE_SUFFIX} --date_range ${DATE_RANGE} --latitude_range ${LATITUDES} --longitude_range ${LONGITUDES}"


python met_extraction_script.py ${ARGUMENTS_1}


