#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import wget
from pathlib import Path
import pyreadr
import datetime
import pandas as pd
import numpy as np

from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn import preprocessing

import geopy.distance as distance
import argparse


if __name__ == '__main__':
    global VERBOSE
    DEFAULT_METADATA_FILE = "AURN_metadata.RData"
    DEFAULT_VERBOSE = 0

    # read arguments from the command line
    parser = argparse.ArgumentParser(description="*** A script for automated downloading of AURN data for a given date range. ***")
    parser.add_argument("--meta_data_url", "-m", help="url of the AURN metadata")
    parser.add_argument("--meta_data_filename", "-f", help="filename of the AURN metadata in RData format (.RData). \
                                                           Default: {}".format(DEFAULT_METADATA_FILE))
    parser.add_argument("--emep_filename","-e", default=None, help="filename of the emep file in CSV format (.csv)")
    parser.add_argument("--years", "-y", metavar='Y', type=int, nargs='+', help="the years to be processed. Must be \
        in (and defaults to) {}".format('[' + ", ".join([str(int) for int in AVAILABLE_YEARS]) + ']'))
    parser.add_argument("--min_years", "-n", type=int, help="minimum number of years of data that a site must have")
    parser.add_argument("--useful_num_years", "-u", type=int, help="minimum number of years of data for any site that \
        we are going to use as a reference site later")
    parser.add_argument("--sites", "-s", metavar='S', dest="sites", type=str, nargs='+', help="the measurement sites \
        to be processed. Default is to process all available AURN sites.")

    parser.add_argument("--save_to_csv",dest="save_to_csv",action='store_true',help="save output into CSV format (default).")
    parser.add_argument("--no_save_to_csv",dest="save_to_csv",action='store_false',help="don't save output to CSV format")
    parser.set_defaults(save_to_csv=True)

    parser.add_argument("--load_from_csv",dest="load_from_csv",action='store_true',help="load input from CSV file.")
    parser.add_argument("--no_load_from_csv",dest="load_from_csv",action='store_false',help="don't load input from csv file (default).")
    parser.set_defaults(load_from_csv=False)

    parser.add_argument("--impute_values",dest="impute_values",action='store_true',help="impute missing values (default).")
    parser.add_argument("--no_impute_values",dest="impute_values",action='store_false',help="don't impute missing values.")
    parser.set_defaults(impute_values=True)

    # Log verbose-ness
    parser.add_argument("--verbose", "-v", type=int,
            help="Level of output for debugging (Default: {} (0 = no verbose output))".format(str(DEFAULT_VERBOSE)))

    # read arguments from the command line
    args = parser.parse_args()

    if args.meta_data_url:
        meta_data_url = args.meta_data_url
    else:
        print('No meta-data_url given, so will look in local folder for meta_data_filename')
        meta_data_url = None

    if args.meta_data_filename:
        meta_data_filename = Path(args.meta_data_filename)
    else:
        print('No meta_data_filename provided, so using default:', DEFAULT_METADATA_FILE)
        meta_data_filename = Path('AURN_metadata.RData')

    if args.emep_filename:
        emep_filename = Path(args.emep_filename)
        if not emep_filename.is_file():
            print('{} does not exist, so not using emep data'.format(emep_filename))
            emep_filename = None
    else:
        print('No emep_filename provided, so not using emep data')
        emep_filename = None

    if args.years:
        years = args.years
        print('Years selected:', years)
    else:
        print('No years provided, so using default: ', '[' + ", ".join([str(int) for int in AVAILABLE_YEARS]) + ']')
        years = AVAILABLE_YEARS

    if args.min_years:
        min_years = args.min_years
        print('Min years (minimum number of years of data that a site must have):', min_years)
    else:
        print('No min_years provided, so using default: 0.4 * number of years: 0.4*{}={}'
              .format(len(years), 0.4*len(years)))
        min_years = 0.4*len(years)

    if args.useful_num_years:
        useful_num_years = max(args.useful_num_years,min_years)
        print('Useful number of years (minimum number of years of data for any site that we are going to use as a \
            reference site later; this cannot be less than min_years):', useful_num_years)
    else:
        print('No useful_num_years provided, so using default: 0.8 * number of years: 0.8*{}={}'
              .format(len(years), 0.8*len(years)))
        useful_num_years = max(0.8*len(years),min_years)

    if args.sites:
        site_list = args.sites
    else:
        print('No sites provided, so using all available sites in metadata file')
        site_list = None


    # process control flags
    if args.save_to_csv:
        save_to_csv = args.save_to_csv
    else:
        print('No save_to_csv provided, so using default: True')
        save_to_csv = True

    if args.load_from_csv:
        load_from_csv = args.load_from_csv
    else:
        print('No load_from_csv provided, so using default: False')
        load_from_csv = False

    if args.impute_values:
        impute_values = args.impute_values
    else:
        print('No impute_values provided, so using default: True')
        impute_values = True

    if args.verbose:
        VERBOSE = max(args.verbose, 0)
        print('verbose: ', VERBOSE)
    else:
        print('No verbose flag provided, so using default: {}'.format(str(DEFAULT_VERBOSE)))
        VERBOSE = DEFAULT_VERBOSE

    # Does the metadatafile exist?
    if meta_data_filename.is_file():
        print("Meta data file already exists in this directory, will use this")
    else:
        print("Downloading Meta data file")
        wget.download(meta_data_url)


    # Read the RData file into a Pandas dataframe
    metadata = pyreadr.read_r(meta_data_filename.name)


    # If a single year is passed then convert to a list with a single value
    if type(years) is int:
        years = [years]
    print('Years:', years)
    current_year = datetime.datetime.now()


    base_path = Path("AURN_data_download")
    if (base_path.is_dir() is False):
        base_path.mkdir()
    # change this if we need more organised data storage later
    data_path = base_path

    # get list of sites to processextract_site_data
    if not site_list:
        site_list = metadata['AURN_metadata']['site_id'].unique()
    else:
        # Todo - Test that sites provided are correct
        pass
    print('Site list', site_list)

    # create the station location dataset
    stations = metadata['AURN_metadata'][['site_id','latitude', 'longitude','site_name']].drop_duplicates()
    stations = stations.rename(columns={"site_id":"SiteID","latitude":"Latitude","longitude":"Longitude"})
    stations = stations.set_index('SiteID')


    # create a dataframe with the hourly dataset for all stations
    hourly_dataframe = extract_site_data(site_list, metadata, years, data_path, save_to_csv)
    hourly_dataframe = hourly_dataframe.rename(columns={'siteID':'SiteID'})

    # apply some filtering of negative and zero values
    spc_list = ['O3','PM10','PM2.5','NO2','NOXasNO2','SO2']
    for spc in spc_list:
        print('filtering {}:'.format(spc))
        try:
            print('\t{} has {} positive values'.format(spc,len(hourly_dataframe.loc[hourly_dataframe[spc]>0.0])))
            print('\t{} has {} NaNs'.format(spc,len(hourly_dataframe.loc[hourly_dataframe[spc].isna()])))
            print('\t{} has {} negative or zero values that will be replaced with NaNs'.format(spc,len(hourly_dataframe.loc[hourly_dataframe[spc]<=0.0])))
            hourly_dataframe.loc[hourly_dataframe[spc]<=0.0, spc] = np.nan
        except:
            print('\t{} has  no values'.format(spc))


    # load the EMEP model data, or create an empty dataframe (required for logic checks in the workflow)
    if emep_filename:
        print('reading emep file')
        emep_dataframe = pd.read_csv(emep_filename)
        emep_dataframe = emep_dataframe.rename(columns={'NOx':'NOXasNO2'})
    else:
        emep_dataframe = pd.DataFrame()

    # pull out the daily mean and max values for the site list
    # postprocessing the data set, to get daily data
    daily_dataframe = postprocess_organisation(hourly_dataframe, emep_dataframe, stations, site_list, impute_values,
        useful_num_years, min_years)


    # sort the data
    daily_dataframe = daily_dataframe.sort_index()

    # write this dataset to file
    daily_dataframe.to_csv(data_path.joinpath('pollution_daily_data_{}-{}.csv'.format(years[0],years[-1])),index=True,header=True,float_format='%.2f')




















