# AURN data download tools

This directory contains the scripts for downloading AURN data for 
a given list of sites and years.

Output data is a csv file with the following columns:
- timestamp
- site_id
- NO2
- SO2
- NOXasNO2
- PM2.5
- PM10
- O3

## Data Extraction

The AURN data download and extraction script is `AURN_download.py`. The required flags
for running this script can be obtained by passing the `--help` flag.

The bash script `RUN_AURN_download.sh` contains settings for two setups, one is to download
all available AURN data for the years 2016-2019 ('FULL'), and the other is to download a 
small subset of stations for the years 2016-2017 ('REDUCED'). These both will extract the
AURN data from the RData files, and save this as a csv data file in the output directory.