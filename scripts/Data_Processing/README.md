# Meteorological and Pollution processing tools

This directory contains the scripts for processing the meteorological and pollution
datasets, to produce the final daily mean and maximum datasets. This processing can
include the imputation of missing data, with control flags to determine the minimum
years of data needed for measurement sites to be included in the daily dataset.

The scripts in the `AURN_Data_Download` and `MEDMI_Data_Download` directories will
need to be run before these are (and the data copied across from the MEDMI server).

The meteorological scripts are:
- `Met_Temperature_RH_Pressure_data_processing.py`
- `RUN_met_script.sh`
The meteorological processing script requires the `station_data_clean.csv` file, which
is used for the station location information (latitude and longitude), and taken from
the MEDMI database. This station data was complete in 2020, but if new stations are
added to the network after this date, it may need updating too.

The AURN pollution scripts are:
- `AURN_processing.py`
- `RUN_aurn_script.sh`
The AURN processing script requires the `AURN_metadata.RData` file, which is downloaded
from the AURN website.

## Common Settings

The settings common to scripts are:
- `--date_range DATE_RANGE`: the start and end dates for the period of interest
- `--min_years MIN_YEARS`: the minimum number of years of data a site must have to be included
      in the final dataset
- `--min_years_ref MIN_YEARS_REF`: the minimum number of years of data a site must have to
      be used as a reference site for imputation processes
- `--impute_values`: use imputation for missing values
- `--no_impute_values`: do not use imputation to fill missing values
- `--verbose VERBOSE`: verbosity of log messages (0=minimal)

## Meteorological Settings

Meteorological specific settings