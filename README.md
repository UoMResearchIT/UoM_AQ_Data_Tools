# UoM_AQ_Data_Tools
This repository contains tools for obtaining and processing UK air quality data.

## Obtaining DEFRA AURN data

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
 - create -n webscraper selenium

The script has been developed on OSX 10.14.6, using chrome 80.0.x and selenium 3.141.0
It was developed in January 2020 - if changes are made to the DEFRA AURN website after
this then the script may need modifying.

To use the script:
- edit the settings at the end of the script (in the labelled settings section)
- python defra_website_automation.py

Do not attempt to do anything else on your computer while it is running - otherwise you
risk disrupting the automated process.

## Processing DEFRA AURN data


## Obtaining Met & Pollen data


## Processing Met & Pollen data


## Combining Met, Pollen, and AURN data
