#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("..")
from pathlib import Path
import argparse


from environmental_data_modules import AurnPostProcessor


if __name__ == '__main__':

    # read arguments from the command line
    parser = argparse.ArgumentParser(description="*** A script for automated downloading of AURN data for a given date range. ***")
    parser.add_argument("--meta_data_url", "-m", help="url of the AURN metadata")
    parser.add_argument("--meta_data_filename", "-f", help="filename of the AURN metadata in RData format (.RData). \
                                                           Default: {}".format(AurnPostProcessor.DEFAULT_METADATA_FILE))
    parser.add_argument("--emep_filename","-e", default=None, help="filename of the emep file in CSV format (.csv)")
    parser.add_argument("--years", "-y", metavar='Y', type=int, nargs='+', help="the years to be processed. Must be \
        in (and defaults to) {}".format('[' + ", ".join([str(int) for int in AurnPostProcessor.AVAILABLE_YEARS]) + ']'))
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
            help="Level of output for debugging (Default: {} (0 = no verbose output))".format(str(AurnPostProcessor.DEFAULT_VERBOSE)))

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
        print('No useful_num_years provided, so using default: max of min_years and 0.8 * number of years: 0.8*{}={})'
              .format(min_years, 0.8*len(years)))
        useful_num_years = max(0.8*len(years), min_years)

    if args.sites:
        site_list = args.sites
        print('Site list: {}'.format(site_list))
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
        print('No verbose flag provided, so using default: {}'.format(str(AurnPostProcessor.DEFAULT_VERBOSE)))
        VERBOSE = AurnPostProcessor.DEFAULT_VERBOSE

    # If a single year is passed then convert to a list with a single value
    if type(years) is int:
        years = [years]
    print('Years:', years)

    processor = AurnPostProcessor(meta_data_url=meta_data_url,
                                  out_dir=AurnPostProcessor.DEFAULT_OUT_DIR, site_list=site_list)
    processor.process(file_in=meta_data_filename, outfile_suffix='_test1', years=years,
                      exclude_site_list=AurnPostProcessor.DEFAULT_EXCLUDE_STATION_LIST,
                      emep_filename=emep_filename,
                      useful_num_years=useful_num_years,
                      save_to_csv=AurnPostProcessor.DEFAULT_SAVE_TO_CSV)


