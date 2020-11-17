#!/bin/bash

OUTDIR_PREFIX='dataxxx_'
OUTFILE_SUFFIX='2017-01--2017-06'
DATE_RANGE='2017-1-1_0 2017-06-30_23'
LATITUDES='53 55'
LONGITUDES='-5 -3'
VERBOSE=1
MEASUREMENTS='temperature rain wind pollen'


# EXAMPLE 1
ARGUMENTS_1="--outdir_prefix ${OUTDIR_PREFIX} --outfile_suffix ${OUTFILE_SUFFIX} --date_range ${DATE_RANGE} --latitude_range ${LATITUDES} --longitude_range ${LONGITUDES}"

# TEST WITH DOUGS
OUTDIR_TEST='Anns_for_Doug_tests'
OUTFILE_SUFFIX_TEST='June2017'
DATE_RANGE_TEST='2017-06-01_0' '2017-06-30_23'

ARGUMENTS_TEST="-o ${OUTDIR_TEST} -s ${OUTFILE_SUFFIX_TEST} -d ${DATE_RANGE_TEST} -m ${MEASUREMENTS} -x -v ${VERBOSE}"

python met_extraction_script.py ${ARGUMENTS_TEST}
