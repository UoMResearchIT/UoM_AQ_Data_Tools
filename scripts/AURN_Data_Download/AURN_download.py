#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Program for downloading the AURN measurement datasets
for a given time period and set of stations. These
are obtained from the RData formatted datasets.

Usage can be obtained using the --help flag.


@author:    Douglas Lowe and Ann Gledson
            Research IT, University of Manchester

"""

import sys
sys.path.append("../..")
import argparse


from environmental_data_modules import AurnExtractor, AurnPostProcessor


if __name__ == '__main__':

    # read arguments from the command line
    parser = argparse.ArgumentParser(description="Downloading script for AURN data for a given date range and sites.")
    parser.add_argument("--metadata_url", "-m",
                        help="url of the AURN metadata. Default: {}".format(AurnPostProcessor.DEFAULT_METADATA_URL))
    parser.add_argument("--metadata_filename", "-f", type=str,
                        help="filename of the AURN metadata in RData format (.RData). " \
                             "Default: {}".format(AurnPostProcessor.DEFAULT_METADATA_FILE))
    parser.add_argument("--years", "-y", metavar='Y', type=int, nargs='+', help="the years to be processed. Must be \
        in (and defaults to) {}".format('[' + ", ".join([str(int) for int in AurnExtractor.get_available_years()]) + ']'))
    parser.add_argument("--sites", "-i", metavar='S', dest="sites", type=str, nargs='+', help="the measurement sites \
        to be processed. Default is to process all available AURN sites.")
    parser.add_argument("--save_to_csv",dest="save_to_csv",action='store_true',help="save output into CSV format (default).")
    parser.add_argument("--no_save_to_csv",dest="save_to_csv",action='store_false',help="don't save output to CSV format")
    parser.set_defaults(save_to_csv=True)

    # output directory/file names
    parser.add_argument("--outdir_name", "-o", dest="outdir_name", type=str,
                        help="output directory name. Default: {}".format(AurnExtractor.DEFAULT_OUT_DIR))
    parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
                        help="suffix to be appended to output file name. \
                              Base name is: {}. Default suffix is: '{}'".format(
                            AurnExtractor.BASE_FILE_OUT,AurnExtractor.DEFAULT_OUT_FILE_SUFFIX))
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

    if args.years:
        years = args.years
        print('Years selected:', years)
    else:
        print('No years provided, so using default: [{}]'.format(", ".join([str(int) for int in AurnExtractor.get_available_years()])))
        years = AurnExtractor.get_available_years()

    if args.sites:
        site_list = args.sites
        print('Site list: {}'.format(site_list))
    else:
        print('No sites provided, so using all available sites in metadata file')
        site_list = None

    print('Save to csv: {}'.format(args.save_to_csv))

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
