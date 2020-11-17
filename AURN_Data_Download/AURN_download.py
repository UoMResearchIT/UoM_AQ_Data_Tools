#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Sep 22 09:16:40 2020

Script for automated downloading of AURN data for a given date range.
Summary of how many sites we will be imputing for data for each species,
and how many useful site there are for each.

NOXasNO2 has 157 req sites and 118 useful sites
O3 has 75 req sites and 66 useful sites
NO2 has 157 req sites and 118 useful sites
PM2.5 has 81 req sites and 55 useful sites
PM10 has 77 req sites and 51 useful sites
SO2 has 27 req sites and 19 useful sites

Todo: Doug following needs updating
We will only impute data for sites with > 1 year of data, and only use
sites with >3.5 years of data for the imputation inputs

Full dataset information on NaN's and negative/zero numbers:
O3 has 2507061 positive values
O3 has 3139956 NaNs
O3 has 735 negative or zero values that will be replaced with NaNs
PM10 has 2287338 positive values
PM10 has 3352010 NaNs
PM10 has 8404 negative or zero values that will be replaced with NaNs
PM2.5 has 2336213 positive values
PM2.5 has 3280461 NaNs
PM2.5 has 31078 negative or zero values that will be replaced with NaNs
NO2 has 4936512 positive values
NO2 has 707558 NaNs
NO2 has 3682 negative or zero values that will be replaced with NaNs
NOXasNO2 has 4939775 positive values
NOXasNO2 has 707425 NaNs
NOXasNO2 has 552 negative or zero values that will be replaced with NaNs
SO2 has 837779 positive values
SO2 has 4806980 NaNs
SO2 has 2993 negative or zero values that will be replaced with NaNs
(note, the NaN count will include sites that we will not be imputing with data)

@author:    Douglas Lowe and Ann Gledson
            Research IT, University of Manchester

"""

import sys
sys.path.append("..")
import argparse


from environmental_data_modules import AurnExtractor, AurnPostProcessor


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
        in (and defaults to) {}".format('[' + ", ".join([str(int) for int in AurnExtractor.get_available_years()]) + ']'))
    parser.add_argument("--min_years", "-n", type=int, help="minimum number of years of data that a site must have")
    parser.add_argument("--min_years_ref", "-u", type=int, help="minimum number of years of data for any site that \
        we are going to use as a reference site later. (this cannot be less than min_years)")
    parser.add_argument("--sites", "-i", metavar='S', dest="sites", type=str, nargs='+', help="the measurement sites \
        to be processed. Default is to process all available AURN sites.")
    parser.add_argument("--save_to_csv",dest="save_to_csv",action='store_true',help="save output into CSV format (default).")
    parser.add_argument("--no_save_to_csv",dest="save_to_csv",action='store_false',help="don't save output to CSV format")
    parser.set_defaults(save_to_csv=True)
    parser.add_argument("--impute_values",dest="impute_values",action='store_true',help="impute missing values (default).")
    parser.add_argument("--no_impute_values",dest="impute_values",action='store_false',help="don't impute missing values.")
    parser.set_defaults(impute_values=True)

    # output directory/file names
    parser.add_argument("--outdir_name", "-o", dest="outdir_name", type=str,
                        help="output directory name. Default: {}".format(AurnExtractor.DEFAULT_OUT_DIR))
    parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
                        help="suffix to be appended to output file name. Default: {}".format(
                            AurnExtractor.DEFAULT_OUT_FILE_SUFFIX))
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
        print("No metadata_filename provided, so using default: '{}".format(AurnPostProcessor.DEFAULT_METADATA_FILE))
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
        print('No years provided, so using default: [{}]'.format(", ".join([str(int) for int in AVAILABLE_YEARS])))
        years = AVAILABLE_YEARS

    if args.min_years:
        min_years = args.min_years
        print('Min years (minimum number of years of data that a site must have):', min_years)
    else:
        print('No min_years provided, so using default: 0.4 * number of years: 0.4*{}={}'
              .format(len(years), 0.4*len(years)))
        min_years = 0.4*len(years)

    if args.min_years_ref:
        min_years_ref = max(args.min_years_ref, min_years)
        print('Minimum years ref:', min_years_ref)
    else:
        min_years_ref = max(min_years, 0.8 * len(years))
        print('No min_years_ref provided, so using default: max(min_years, 0.8*number of years): \
         Max({}, 0.8*{}={})'.format(min_years, len(years), max(min_years, 0.8*len(years))))

    if args.sites:
        site_list = args.sites
        print('Site list: {}'.format(site_list))
    else:
        print('No sites provided, so using all available sites in metadata file')
        site_list = None

    print('Save to csv: {}'.format(args.save_to_csv))
    print('Impute values: {}'.format(args.impute_values))

    if args.outdir_name:
        outdir_name = args.outdir_name
        print('Using outdir_name: {}'.format(outdir_name))
    else:
        print('No outdir_name given, so will use default: {}'.format(AurnExtractor.DEFAULT_OUT_DIR))
        outdir_name = AurnExtractor.DEFAULT_OUT_DIR

    if args.outfile_suffix:
        outfile_suffix = args.outfile_suffix
        print('Using outfile_suffix: {}'.format(outfile_suffix))
    else:
        print('No outfile_suffix provided, so using default: {}'.format(str(AurnExtractor.DEFAULT_OUT_FILE_SUFFIX)))
        outfile_suffix = AurnExtractor.DEFAULT_OUT_FILE_SUFFIX

    if args.verbose:
        verbose = max(args.verbose, 0)
        print('verbose: ', verbose)
    else:
        print('No verbose flag provided, so using default: {}'.format(str(AurnPostProcessor.DEFAULT_VERBOSE)))
        verbose = AurnPostProcessor.DEFAULT_VERBOSE

    extractor = AurnExtractor(metadata_filename=metadata_filename,
                              metadata_url=metadata_url,
                              out_dir=outdir_name, verbose=verbose)
    extractor.extract_data(
                           years=years,
                           site_list=site_list,
                           save_to_csv=args.save_to_csv,
                           outfile_suffix=outfile_suffix)

    processor = AurnPostProcessor(metadata_filename=metadata_filename,
                                  metadata_url=metadata_url,
                                  out_dir=outdir_name, verbose=verbose)
    processor.process(  in_file=extractor.file_out,
                        outfile_suffix= outfile_suffix,
                        site_list=site_list,
                        emep_filename=emep_filename,
                        min_years=min_years,
                        min_years_reference=min_years_ref,
                        impute_data=args.impute_values,
                        save_to_csv=args.save_to_csv)
