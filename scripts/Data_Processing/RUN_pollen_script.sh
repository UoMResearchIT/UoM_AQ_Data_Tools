#!/bin/bash

## Fixed arguments. Adjust the verbosity (0=none) as required.
VERBOSE=2

FIXED_ARGS="--verbose ${VERBOSE}"


## Scenario settings
DATA_DIR="../MEDMI_Data_Download/full_data"
OUTFILE_SUFFIX="2016-2019"


SCENARIO_ARGS="--outfile_suffix ${OUTFILE_SUFFIX} --data_dir ${DATA_DIR}"


python Pollen_data_processing.py ${FIXED_ARGS} ${SCENARIO_ARGS}