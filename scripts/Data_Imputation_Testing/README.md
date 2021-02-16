# Meteorological and Pollution imputation testing tools

This directory contains the scripts for testing the imputation of the meteorological 
and pollution datasets. This is not part of the main workflow for producing the final
dataset, but can be used to investigate the implications of the assumptions made in
that workflow.

The scripts will remove data from a given fraction of the time period requested, for 
the sites to be used for testing. These sites must have enough data to satisfy the
reference site requirement, a list of suitable sites for your dataset and settings 
can be obtained using the `--check_sites` flag. The position of the data to be removed
for each site tested can be specified, either as contiguous blocks of data at the start,
middle, or end of the timeseries, or as randomly distributed losses within the whole
timeseries. The sites with data removed will be tested individually, and may be used
as reference sites (without their data removed) for other sites being tested.

The scripts in the `AURN_Data_Download` and `MEDMI_Data_Download` directories will
need to be run before these are (and the data copied across from the MEDMI server).

The meteorological scripts are:
- `Met_Temperature_RH_Pressure_impute_testing.py`
- `RUN_met_impute_script.sh`
- `met_stat_plotting.py`
The meteorological processing script requires the `station_data_clean.csv` file, which
is used for the station location information (latitude and longitude), and taken from
the MEDMI database. This station data was complete in 2020, but if new stations are
added to the network after this date, it may need updating too.

The AURN pollution scripts are:
- `AURN_impute_testing.py`
- `RUN_aurn_impute_script.sh`
- `aurn_stat_plotting.py`
The AURN processing script requires the `AURN_metadata.RData` file, which is downloaded
from the AURN website.

## Data Analysis

The statistical data for each test (as a csv file), as well as plots of the comparisons 
of the hourly and daily data (as pdf files), will be stored within local `met_stats` 
or `aurn_stats` directories (unless the output paths are changed). The aggregated stats can
be plotted using the `met_stat_plotting.py` or `aurn_stat_plotting.py` scripts. These 
are examples for the data paper, and will need adapting to fit the analysis simulations you
are interested in.


## Common Settings

The settings common to scripts are:
- `--check_sites`: list sites suitable for testing imputation (based on settings below), 
      then exits script.
- `--data_lost DATA_LOST`: fraction of the total time period to remove data for. 
- `--data_lost_position POSITION`: position of data to be removed (`start`, `middle `, `end`, or `random`) 
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
