# MEDMI data extraction tools

This directory contains the scripts for extracting the MEDMI meteorological measurement
data for given geographic and time ranges. 

This must be run on the MEDMI server. Details of how to do this can be obtained from 
the MED-MI website: [https://www.data-mashup.org.uk/contact-us/](https://www.data-mashup.org.uk/contact-us/)
and  following the instructions for connecting to the server using an SSH client.
Alternatively, email `health@metoffice.gov.uk` to request access.

Output data is a set of csv files, one for each of the requested datasets, each with a 
header listing the dataset extracted and date range, with the following two fixed columns:
timestamp and site_id. The other columns are determined by the datasets.

Available dataset:
- temperature
- pressure
- rel_hum (relative humidity)
- dewpoint (dewpoint temperature)
- rain
- wind
- pollen (selects all the following pollen subsets)
  - salix
  - ambrosia
  - fraxinus
  - ulmus
  - artemisia
  - platanus
  - urtica
  - betula
  - poaceae
  - alnus
  - corylus
  - quercus




## Data Extraction

The data extraction script is `met_extraction_script.py`. This requires the following flags to
be passed:
- `--measurements [LIST]`: is a list of measurements to extract, from the above datasets
  - default is to extract all datasets listed above
- `--extra_measurements`: flag to switch on the extraction of extra met datasets for the requested 
  datasets. This only affects the following datasets (and will extract them all):
  - `temperature`
  - `pressure`
  - `rel_hum` (relative humidity)
  - `dewpoint` (dewpoint temperature)  
- `--outdir_name OUTDIR_NAME`: this is the output directory to use for the final datafiles.
- `--outfile_suffix OUTFILE_SUFFIX`: this is a suffix to append to the output file name.
- `--date_range DATE_RANGE`: this should be a pair of start and end dates, following this format:
  - `YYYY-MM-DD_H[H]`
- `--latitude_range LAT_RANGE`: this should be a pair of values giving the minimum and maximum latitudes.
- `--longitude_range LONG_RANGE`: this should be a pair of values giving the minimum and maximum longitudes.
- `--verbose VERBOSE`: a numerical value indicating the verbosity of the running logs (default is 0 - minimal)
  
The `bash_medmi_extraction.sh` script has two examples for extracting the required measurements.
The reduced scenario is a limited area and period subset of the full scenario (which covers the
whole UK, and the years 2016-2019). The extracted data will be stored as a set of files in
the output directory, following the naming convention: 
`Met_extracted_{DATASET}-{EXTRA_DATASETS}_{OUTFILE_SUFFIX}.csv`


