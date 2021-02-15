# EMEP data extraction tools

This directory contains the scripts for extracting hourly EMEP model data for
a given set of locations.

Required data:
- EMEP output files
- A single WRF output file, for the same domain
- station metadata file (either `RData` or `csv` format)
  - this can be obtained by running the AURN download scripts

Output data is a csv file with the following columns:
- timestamp
- site_id
- NO2
- SO2
- NOx
- PM2.5
- PM10
- O3

## Data Extraction

The data extraction script is `emep_extract.py`. This requires the following flags to
be passed:
- `--sensor_file SENSOR_FILE`: this is the path to the sensor/station metadata file
- `--sensor_file_type [CSV|RDATA]`: specify if the file is in `csv` or `RData` format
  - for `csv` format the columns expected are `long`, `lat`, and `sensor_name`
  - for `RData` format the AURN metadata format is expected, containing `longitude`, 
    `latitude`, and `site_id`.
- `--wrf_file WRF_FILE`: this is the path to the WRF output file (for location information)
  - this only has to be the same grid as the EMEP data, it does not have to be for the same
    time period.
- `--emep_file EMEP_FILE`: this is the path to the EMEP file to be processed.
- `--out_file OUT_FILE`: this is the path and name to use for the output csv file.
- flag to indicate the type of output expected:
  - `--daily_means`: output daily mean/max data
  - `--no_daily_means`: output the raw data (probably hourly)
  
The `batch_hourly_emep_extraction.sh` script is an example SGE batch array script used for 
extracting data from 24 EMEP files. This outputs a separate csv file for each period, in
the directory `extracted_model_data_hourly` (following the example script).

## Data Combination 

The individual period data files can then be combined using the script `combine_data.sh`
in the directory `extracted_model_data_hourly`. This script removes the last day from each 
period, to avoid overlaps with the next period - check your data files to determine if this 
is necessary or not.

