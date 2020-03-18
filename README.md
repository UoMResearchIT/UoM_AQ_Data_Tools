# UoM_AQ_Data_Tools
This repository contains tools for obtaining and processing UK air quality data.

## Obtaining DEFRA AURN Measurement Data

`defra_website_automation.py`

This script will automate the filling in and submitting of the DEFRA AURN website. 
This generates requests for data, which will be emailed (to the email address you 
provide in the script) once they have been extracted.

To run the script you need to install Selenium and Chrome, there are instructions to do
this on this webpage: https://selenium-python.readthedocs.io/installation.html

In summary (and using conda, not pure pip), install:
- chrome
- chrome drivers for selenium (https://sites.google.com/a/chromium.org/chromedriver/downloads)
- selenium, via conda, creating a virtual environment for this:
 - conda create -n webscraper selenium

The script has been developed on OSX 10.14.6, using chrome 80.0.x and selenium 3.141.0
It was developed in February 2020 - if changes are made to the DEFRA AURN website after
this then the script may need modifying.

To use the script:
- edit the settings at the end of the script (in the labelled settings section)
- python defra_website_automation.py

Do not attempt to do anything else on your computer while it is running - otherwise you
risk disrupting the automated process.

## Obtaining DEFRA AURN Station Location Information

Note: this step only needs to be followed if more AURN stations have been added since the
time of writing (March 2020). Otherwise the included csv data file can be used. 

To download a list of all DEFRA AURN stations use this website:
https://uk-air.defra.gov.uk/networks/find-sites?view=advanced
Select the "Include closed monitoring sites in search" radio button, leave all other 
options as the defaults (giving as wide a search as possible), and press the green 
"search network" button.

The station information contains latitude and longitude, but no address. To get the address
information use this script `aurn_address_postcode_calculator.py`. This will require you
to create a Google API key for the Geocoding service - you can do this here: 
https://console.developers.google.com (at the time of writing, March 2020, this is free
for reasonable usage in the first year).

This script creates the file `aurn_measurement_sites_addresses_postcodes.csv`. A copy of
this file, containing information for 272 AURN stations is included in this repository.


## Processing DEFRA AURN data

The raw data files obtained from the DEFRA website can be processed using the script
`aurn_files_combined_and_process.py`. This invoked as:
`python aurn_files_combine_and_process.py -i [input directory] -o [output directory] -s [AURN station file]`
All input flags are optional, if they are not included then default values are used
(the current working directory for input / output directory, and the file 
`aurn_measurement_sites_addresses_postcodes.csv` for the AURN station file).

A single csv file is generated, following this naming convention:
`AURN_AQ_data_[start date]_[end date].csv`
The start and end dates will be for the data contained within that file. 

Note that the date range of the final file may shorter than the date ranges that you 
have downloaded, as this script drops all non-validated data points (and data validation 
can take several months to perform). If you wish to include all data (validated or otherwise) 
then you can comment out this line of code in the function `extract_aurn_values_data`: 
`value_data[bool_mask] = np.nan` 


## Obtaining Met & Pollen data


## Processing Met & Pollen data


## Combining Met, Pollen, and AURN data
