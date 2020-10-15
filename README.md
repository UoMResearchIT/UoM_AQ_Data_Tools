# UoM AQ Data Tools

This repository contains tools for obtaining and processing UK air quality data.

The sections below are:
 - Requirements
 - Extracting EMEP Model Air Quality Data
 - Obtaining and Processing AURN Air Quality Data
 - Obtaining MEDMI Meteorological and Pollen Data
 - Processing MEDMI Data
 - Combining Datasets


## Requirements

The processing scripts in this repository are written in python, tested on unix and OSX
systems.

The EMEP (and WRF) models are written in fortran - follow these references for compiling 
theses, preparing model inputs, and performing the simulations.

The MEDMI dataset are accessed using the python2 installation on their system, no more
packages require installing to run the scripts for this.

The processing scripts for extracting the EMEP data are written in python3. To
install the packages needed for these (using conda and pip) use this script:
`conda env create -f env_emep.yml`

The processing scripts for obtaining the AURN dataset, and processing all datasets, are
written in python3. To install the packages needed for these, use this script: 
`conda env create -f env_aurn_medmi.yml`


## Extracting EMEP Model Air Quality Data

EMEP simulations should be carried out using the tools available in the X repository.

Scripts for extracting the required data are given in the `EMEP_Data_Processing` directory.

These will extract the hourly data, and daily mean / max values, for locations taken
from the AURN metadata files. Bash scripts are provided for running these tools.
 

## Obtaining and Processing AURN Air Quality Data

(All scripts for processing DEFRA AURN data are in the `AURN_DATA_Download` directory)

`AURN_download.py`

This script accesses the AURN datasets that are provided in RData format on the UK-Air
DEFRA website.

Run using python and use the -help command to get the full set of parameters/options.

This script creates a `AURN_data_download` directory, in which the available RData format
files for all AURN sites will be downloaded and stored, and then will create a combined
data file called `pollution_daily_data_[start year]-[end year].csv` in the same directory.
This data file contains the combined daily mean, max, and data count for all AURN sites.


## Obtaining MEDMI Meteorological and Pollen Data

All scripts for obtaining Met and Pollen data are in the `MEDMI_Data_Download` directory.
The scripts are:
- `met_extraction_script.py`
- `pollen_extraction_script.py`
- `pollen_site_met_calculations.py`


To use these scripts, it is neccessary to log onto the MEDMI ssh server. Details of how to do 
this can be obtained from the MED-MI website:
https://www.data-mashup.org.uk/contact-us/
and  following the instructions for:  ***connect to the server using an SSH client***.
Alternatively, email: *health at metoffice dot gov dot uk*

Once logged in, create these directories for the extracted data : `data_met`, `data_rain`, and `data_pollen`
(to change these directory names, in the scripts edit the `file_base` or `fname` strings).

Run using python and use the -help command to get the full set of parameters/options.

Copy the scripts to the ssh server. Then run each in turn (simply `python [script name] [params/options]`).
The `pollen_site_met_calculations.py` script will take longest to run, as this interpolates
met data to the pollen measurement site locations - it is computationally intensive and,
depending on the time period you've chosen, could take several hours to complete.

Once all scripts have finished running you can copy the data back to your local computer.

## Processing MEDMI Data


## Combining Datasets
