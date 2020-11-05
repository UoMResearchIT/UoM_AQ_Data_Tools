#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("..")
import argparse


from environmental_data_modules import AurnPostProcessor


if __name__ == '__main__':

    # read arguments from the command line
    parser = argparse.ArgumentParser(description="*** A script for automated downloading of AURN data for a given date range. ***")
    parser.add_argument("--metadata_url", "-m",
                        help="url of the AURN metadata. Default: {}".format(AurnPostProcessor.DEFAULT_METADATA_URL))
    parser.add_argument("--metadata_filename", "-f", type=str,
                        help="filename of the AURN metadata in RData format (.RData). " \
                             "Default: {}".format(AurnPostProcessor.DEFAULT_METADATA_FILE))
    parser.add_argument("--emep_filename","-e", default=None, help="filename of the emep file in CSV format (.csv)")
    parser.add_argument("--years", "-y", metavar='Y', type=int, nargs='+', help="the years to be processed. Must be \
        in (and defaults to) {}".format('[' + ", ".join([str(int) for int in AurnPostProcessor.AVAILABLE_YEARS]) + ']'))
    parser.add_argument("--min_years", "-n", type=int, help="minimum number of years of data that a site must have")
    parser.add_argument("--useful_num_years", "-u", type=int, help="minimum number of years of data for any site that \
        we are going to use as a reference site later. (this cannot be less than min_years)")
    parser.add_argument("--sites", "-s", metavar='S', dest="sites", type=str, nargs='+', help="the measurement sites \
        to be processed. Default is to process all available AURN sites.")

    parser.add_argument("--save_to_csv",dest="save_to_csv",action='store_true',help="save output into CSV format (default).")
    parser.add_argument("--no_save_to_csv",dest="save_to_csv",action='store_false',help="don't save output to CSV format")
    parser.set_defaults(save_to_csv=True)

    parser.add_argument("--impute_values",dest="impute_values",action='store_true',help="impute missing values (default).")
    parser.add_argument("--no_impute_values",dest="impute_values",action='store_false',help="don't impute missing values.")
    parser.set_defaults(impute_values=True)

    # Log verbose-ness
    parser.add_argument("--verbose", "-v", type=int,
            help="Level of output for debugging (Default: {} (0 = no verbose output))".format(
                AurnPostProcessor.DEFAULT_VERBOSE))

    # read arguments from the command line
    args = parser.parse_args()

    if args.metadata_url:
        metadata_url = args.metadata_url
        print('Metadata url: {}'.format(metadata_url))
    else:
        print('No metadata_url given, so will use {}, if no metadata filename provided.'.format(
            AurnPostProcessor.DEFAULT_METADATA_URL))
        metadata_url = None

    if args.metadata_filename:
        metadata_filename = args.metadata_filename
        print('Metadata filename: {}'.format(metadata_filename))
    else:
        print('No metadata_filename provided, so using default:', AurnPostProcessor.DEFAULT_METADATA_FILE)
        metadata_filename = AurnPostProcessor.DEFAULT_METADATA_FILE

    if args.emep_filename:
        emep_filename = args.emep_filename
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
        print('Useful number of years:', useful_num_years)
    else:
        useful_num_years = max(min_years, 0.8 * len(years))
        print('No useful_num_years provided, so using default: max(min_years, 0.8*number of years): \
         Max({}, 0.8*{}={})'.format(min_years, len(years), max(min_years, 0.8*len(years))))


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
        save_to_csv = True
    print('Save to csv: {}'.format(save_to_csv))

    if args.impute_values:
        impute_values = args.impute_values
    else:
        impute_values = True
    print('Impute values: {}'.format(impute_values))

    if args.verbose:
        verbose = max(args.verbose, 0)
        print('verbose: ', verbose)
    else:
        print('No verbose flag provided, so using default: {}'.format(str(AurnPostProcessor.DEFAULT_VERBOSE)))
        verbose = AurnPostProcessor.DEFAULT_VERBOSE


    processor = AurnPostProcessor(out_dir=AurnPostProcessor.DEFAULT_OUT_DIR, verbose=verbose)
    processor.process(metadata_filename=metadata_filename,
                      metadata_url=metadata_url,
                      outfile_suffix='_test1',
                      years=years,
                      site_list=site_list,
                      emep_filename=emep_filename,
                      useful_num_years=useful_num_years,
                      save_to_csv=AurnPostProcessor.DEFAULT_SAVE_TO_CSV)


