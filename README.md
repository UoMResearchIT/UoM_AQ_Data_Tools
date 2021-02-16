# UoM AQ Data Tools

This repository contains tools for obtaining and processing UK air quality data.

The sections below are:
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Data Acquisition](#data-acquisition)
  - [EMEP Model Air Quality Data](#emep-model-air-quality-data)
  - [AURN Air Quality Data](#aurn-air-quality-data)
  - [MEDMI Meteorological and Pollen Data](#medmi-meteorological-and-pollen-data)
- [Data Processing](#data-processing)
  - [MEDMI Data](#medmi-data)
  - [AURN Data](#aurn-data)
  - [Combining Datasets](#combining-datasets)
- [Testing Imputation Methods](#testing-imputation-methods)

<!-- toc -->

## Repository Structure

Operational scripts are stored within the `scripts` directory. Each of the subdirectories 
within this contains a README file guiding tool usage. The `environmental_data_modules` 
directory contains the python modules used by the tools. The `station_data` directory 
contains station specific metadata which has been collated for the datasets, these have 
been gathered from the MEDMI and AURN sources.

```
.
├── environmental_data_modules
├── scripts
│   ├── AURN_Data_Download
│   ├── Combine_Data
│   ├── Data_Imputation_Testing
│   ├── Data_Processing
│   ├── EMEP_Data_Extraction
│   └── MEDMI_Data_Download
└── station_data
```


## Requirements

The processing scripts in this repository are written in python, tested on unix and OSX
systems.

The EMEP (and WRF) models are written in fortran - follow the references below for compiling 
these, preparing model inputs, and performing the simulations.

The MEDMI dataset are accessed using the python2 installation on their system, no more
packages require installing to run the scripts for this.

We recommend using conda to import the required python libraries. Using standalone pip is
untested. 

The processing scripts for extracting the EMEP data are written in python3. To
install the packages needed for these (using conda and pip) use this script:
`conda env create -f env_emep.yml`.
To activate this environment use `conda activate emep`.

The processing scripts for obtaining the AURN dataset, and processing all datasets, are
written in python3. To install the packages needed for these, use this script: 
`conda env create -f env_aurn_medmi.yml`.
To activate this environment this `conda activate aurn_medmi`.

## Data Acquisition

Three separate sets of scripts are provided for obtaining the meteorological and pollution
data. These are covered below and can be run separately as required.

### MEDMI Meteorological and Pollen Data

MEDMI data is obtained using the tools in the `scripts/MEDMI_Data_Download` directory.
These scripts must be run on the MEDMI server, instructions on gaining access to this are
given in the README file in that directory.

### AURN Air Quality Data

Scripts for the downloading of DEFRA AURN data are in the `scripts/AURN_DATA_Download` directory.

These scripts access the AURN datasets that are provided in RData format on the UK-Air
DEFRA website. Guidance on using the scripts is included in the `README.md` file in that
directory, and bash scripts with example configurations are included too.

### EMEP Model Air Quality Data

EMEP simulations should be carried out using the tools available in this Zenodo repository:
[https://doi.org/10.5281/zenodo.3997300](https://doi.org/10.5281/zenodo.3997300).

Scripts for extracting the required data are given in the `scripts/EMEP_Data_Extraction` 
directory. These require use of the `emep` conda environment.

These will extract the hourly data, and daily mean / max values, for locations taken
from the AURN metadata files. Bash scripts are provided for running these tools. The AURN
metadata file will be needed for location information for the extraction of data from the
EMEP files - this can be obtained by running the AURN download scripts first.

## Data Processing

Scripts for processing both the meteorological and pollution datasets are included in
the `scripts/Data_Processing` directory. These require the use of the `aurn_medmi` conda 
environment. Example bash scripts are provided for running these tools.

The processing of the meteorological and pollution datasets is carried out independently of
each other, so can be run as you wish. If EMEP data is to be used in the imputation of the
pollution data, however, that will need to be processed first.

### Combining Datasets

The data files produced using the scripts above are combined using a pandas dataframe merge,
based on date and site identifier. Scripts illustrating how this has been done are included
in the `scripts/Combine_Data` directory.

## Testing Imputation Methods

Imputation of missing data is carried out for the processing of some of the AURN and MEDMI 
datasets. Tools for carrying out statistical testing of the imputation tools used are 
available in the `scripts/Data_Imputation_Testing` directory, with some guidance on using 
these in the README file in that directory.


## Copyright & Licensing

This software has been developed by the [Research IT](https://research-it.manchester.ac.uk/) 
group at the [University of Manchester](https://www.manchester.ac.uk/) for an 
[Alan Turing Institute](https://www.turing.ac.uk/) project.

(c) 2019-2021 University of Manchester.
Licensed under the GPL-3.0 license, see the file LICENSE for details.
