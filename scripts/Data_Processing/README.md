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

The pollen scripts are:
- `Pollen_data_processing.py`
- `RUN_pollen_script.sh`
These require only the pollen data files copied across from the MEDMI server.

## Common Settings for Meteorological and AURN pollution scripts

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

Meteorological specific settings are:
- `--out_dir OUT_DIR`: output directory name
- `--outfile_suffix OUTFILE_SUFFIX`: suffix to append to output filename
- `--file_in FILE_IN`: input file containing the hourly meteorological data to be processed
- `--stations_filename STATIONS_FILENAME`: file containing the station location information
- `--ref_num_stations REF_NUM_STATIONS`: number of reference stations to be used for imputation
- `--min_temp MIN_TEMP`: minimum temperature for measurements (values lower than this are set to NaN)
- `--exclude_sites X`: list of measurement sites to exclude
- `--print_stats`: print statistics for dataset
- `--no_print_stats`: don't print dataset statistics

## AURN pollution Settings

AURN specific settings are:
- `--metadata_url METADATA_URL`: url of AURN metadata
- `--metadata_filename METADATA_FILENAME`: name of AURN metadata file
- `--emep_filename EMEP_FILENAME`: path and filename of the EMEP data file
- `--sites S`: list of the measurement sites to be processed (default is all)
- `--species S`: list of pollutant species to process (default is all)

## Pollen Settings

Pollen specific settings are:
- `--data_dir DATA_DIR`: data directory for both the original files, and the combined file that is produced
- `--file_template FILE_TEMPLATE`: a string giving the standard name for the pollen data files
- `--outfile OUT_FILE`: name of the combined data file that will be output (optional)
- `--outfile_suffix OUT_SUFFIX`: a string giving the identifier to be added to the file template for the output file (optional)
- `--pollen_list POLLENLIST`: a list of the pollen species to process (default is all)
- `--verbose VERBOSE`: verbosity of log messages (0=minimal)