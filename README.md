# UoM AQ Data Tools

This repository contains tools for obtaining and processing UK air quality data.

The sections below are:
- [Requirements](#heading)
 - Extracting EMEP Model Air Quality Data
 - Obtaining and Processing AURN Air Quality Data
 - Obtaining MEDMI Meteorological and Pollen Data
 - Processing MEDMI Data
 - Combining Datasets

```
.
├── environmental_data_modules
├── scripts
│   ├── AURN_Data_Download
│   ├── Data_Imputation_Testing
│   ├── Data_Processing
│   │   └── inputs
│   ├── EMEP_Data_Extraction
│   └── MEDMI_Data_Download
└── station_data
```

<!-- toc -->

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
`conda env create -f env_emep.yml` \
To activate this environment: `conda activate emep`


The processing scripts for obtaining the AURN dataset, and processing all datasets, are
written in python3. To install the packages needed for these, use this script: 
`conda env create -f env_aurn_medmi.yml` \
To activate this environment: `conda activate aurn_medmi`


## Extracting EMEP Model Air Quality Data

EMEP simulations should be carried out using the tools available in this Zenodo repository:
https://doi.org/10.5281/zenodo.3997300.

Scripts for extracting the required data are given in the `scripts/EMEP_Data_Extraction` directory.

These will extract the hourly data, and daily mean / max values, for locations taken
from the AURN metadata files. Bash scripts are provided for running these tools.
 

## Obtaining and Processing AURN Air Quality Data

(All scripts for processing DEFRA AURN data are in the `scripts/AURN_DATA_Download` directory)

`AURN_download.py`

This script accesses the AURN datasets that are provided in RData format on the UK-Air
DEFRA website.

Run using python and use the -help command to get the full set of parameters/options.

This script creates a `AURN_data_download` directory, in which the available RData format
files for all AURN sites will be downloaded and stored, and then will create a combined
data file called `pollution_daily_data_[start year]-[end year].csv` in the same directory.
This data file contains the combined daily mean, max, and data count for all AURN sites.


## Obtaining MEDMI Meteorological and Pollen Data

(All scripts for obtaining Met and Pollen data are in the `scripts/MEDMI_Data_Download` directory)

`met_extraction_script.py`


To use these script(s), it is neccessary to log onto the MEDMI ssh server. Details of how to do 
this can be obtained from the MED-MI website:
https://www.data-mashup.org.uk/contact-us/
and  following the instructions for:  ***connect to the server using an SSH client***.
Alternatively, email: *health at metoffice dot gov dot uk*

Copy the script to the ssh server. Then run it using python (simply `python [script name] [params/options]`) 
and use the --help command to get the full set of parameters/options.

The scripts will create a default directory 'met_extracted_data' for the extracted data.
(To change the directory name from the default use the "--outdir_name" (-o) parameter.)

Once all scripts have finished running you can copy the data back to your local computer.

#### Using extra data
With certain measurements, extra datasets can be requested and added to the outputs. To do this set the 
--extra_measurements (-x) parameter to True (default is False). 
(Note that not all measurements can be extra datasets of others, for example the stations which measure pollen 
don't measure meteorological parameters (and vice versa).)
Therefore, in this tool, the only datasets that can have extra measurements added are:
- temp
- pressure
- dewpoint
- rel_hum

If --extra_measurements is set to True then all allowed extra datasets added will be the above measurements. 
For example temp would have added: pressure, dewpoint and rel_hum


## Processing AURN and MEDMI Data

Scripts for processing these datasets are included in the `scripts/Data_Processing` directory. These require
the use of the `aurn_medmi` conda environment. Example bash scripts are provided for running these. 


## Combining Datasets

TBD

## Testing Imputation Methods

Imputation of missing data is carried out for the processing of some of the AURN and MEDMI datasets. Tools 
for carrying out statistical testing of the imputation tools used are available in the `scripts/Data_Imputation_Testing`
directory.


## Copyright & Licensing

This software has been developed by the [Research IT](https://research-it.manchester.ac.uk/) group at the [University of Manchester](https://www.manchester.ac.uk/) for an [Alan Turing Institute](https://www.turing.ac.uk/) project.

(c) 2019-2021 University of Manchester.
Licensed under the GPL-3.0 license, see the file LICENSE for details.
