#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Sep 22 09:16:40 2020

Script for testing the imputation of AURN data.

@author:    Douglas Lowe and Ann Gledson
            Research IT, University of Manchester

"""

from pathlib import Path
import sys
sys.path.append("../..")
import argparse


from environmental_data_modules import AurnExtractor, AurnImputationTest


if __name__ == '__main__':

    # read arguments from the command line
    parser = argparse.ArgumentParser(description="Script for testing the imputation of AURN data.")
    parser.add_argument("--metadata_url", "-m",
                        help="url of the AURN metadata. Default: {}".format(AurnImputationTest.DEFAULT_METADATA_URL))
    parser.add_argument("--metadata_filename", "-f", type=str,
                        help="filename of the AURN metadata in RData format (.RData). " \
                             "Default: {}".format(AurnImputationTest.DEFAULT_METADATA_FILE))

    ## Dates
    parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',
                        help="start and end dates. (array - first two values only). \
                            Expected date format: {} \n Default: [{}, {}]".format(
                            AurnImputationTest.INPUT_DATE_FORMAT.replace('%', ''),
                            AurnImputationTest.get_available_start(), AurnImputationTest.get_available_end()))

    parser.add_argument("--emep_filename","-e", default=None, help="filename of the emep file in CSV format (.csv)")
    parser.add_argument("--min_years", "-n", type=float, help="minimum number of years of data that a site must have")
    parser.add_argument("--min_years_ref", "-u", type=float, help="minimum number of years of data for any site that \
        we are going to use as a reference site later. (this cannot be less than min_years)")
    parser.add_argument("--sites", "-i", metavar='S', dest="sites", type=str, nargs='+', help="the measurement sites \
        to be processed. Default is to process all available AURN sites.")
    parser.add_argument("--species", "-j", metavar='S', dest="species", type=str, nargs='+', help="the chemical species \
        to be processed. Default is to process all in list: {}.".format(AurnImputationTest.SPECIES_LIST_EXTRACTED))
    parser.add_argument("--save_to_csv",dest="save_to_csv",action='store_true',help="save output into CSV format (default).")
    parser.add_argument("--no_save_to_csv",dest="save_to_csv",action='store_false',help="don't save output to CSV format")
    parser.add_argument("--data_lost", type=float, help="fraction of total period to remove for imputation tests")
    parser.add_argument("--data_loss_position", type=str, help="position to lose data: start, middle, end, random")
    parser.add_argument("--check_sites",dest="check_sites",action='store_true',help="list sites suitable for \
        tests then exit")
    parser.set_defaults(save_to_csv=True)
    parser.set_defaults(impute_values=True)
    parser.set_defaults(check_sites=False)

    # output directory/file names
    parser.add_argument("--outdir_name", "-o", dest="outdir_name", type=str,
                        help="output directory name. Default: {}".format(AurnExtractor.DEFAULT_OUT_DIR))
    parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
                        help="suffix to be appended to output file name. Default: {}".format(
                            AurnExtractor.DEFAULT_OUT_FILE_SUFFIX))
    parser.add_argument("--statdir_name", dest="statdir_name", type=str,
                        help="stat output directory name. Default: {}".format('./'))

    # Log verbose-ness
    parser.add_argument("--verbose", "-v", type=int,
            help="Level of output for debugging (Default: {} (0 = no verbose output))".format(
                AurnImputationTest.DEFAULT_VERBOSE))

    # read arguments from the command line
    args = parser.parse_args()

    if args.data_lost:
        data_lost = args.data_lost
        print('Removing {} fraction of data for imputation tests'.format(data_lost))
    else:
        data_lost = AurnImputationTest.DEFAULT_DATA_LOST
        print('Default: Removing {} fraction of data for imputation tests'.format(data_lost))

    if args.data_loss_position:
        data_loss_position = args.data_loss_position
        print('Using method "{}" for selecting data to lose'.format(data_loss_position))
    else:
        data_loss_position = AurnImputationTest.DEFAULT_DATA_LOSS_POSITION
        print('Default: Using method "{}" for selecting data to lose'.format(data_loss_position))


    if args.metadata_url:
        metadata_url = args.metadata_url
        print('Metadata url: {}'.format(metadata_url))
    else:
        print('No metadata_url given, so will use {}, if no metadata filename provided.'.format(
            AurnImputationTest.DEFAULT_METADATA_URL))
        metadata_url = None

    if args.metadata_filename:
        metadata_filename = args.metadata_filename
        print('Metadata filename: {}'.format(metadata_filename))
    else:
        print("No metadata_filename provided, so using default: '{}".format(AurnImputationTest.DEFAULT_METADATA_FILE))
        metadata_filename = AurnImputationTest.DEFAULT_METADATA_FILE

    if args.emep_filename:
        emep_filename = args.emep_filename
        if not Path(emep_filename).exists():
            print('{} does not exist, so not using emep data'.format(emep_filename))
            emep_filename = None
    else:
        print('No emep_filename provided, so not using emep data')
        emep_filename = None


    if args.date_range:
        if len(args.date_range) >= 2:
            if len(args.date_range) > 2:
                print('Warning: You have input more than 2 dates, only the first 2 will be used.')
            date_range = args.date_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 dates from input --date_range: {}'.format(str(args.date_range)))
        print('Using date range: [{}]'.format(','.join(date_range)))
    else:
        print('No date_range provided, so using default: [{}]'.format(','.join(AurnImputationTest.get_available_dates())))
        date_range = AurnImputationTest.get_available_dates()


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

    if args.species:
        species_list = args.species
        print('Species list: {}'.format(species_list))
    else:
        print('No species provided, so the default list of species in the processor')
        species_list = AurnImputationTest.SPECIES_LIST_EXTRACTED

    if args.outdir_name:
        outdir_name = args.outdir_name
        print('Using outdir_name: {}'.format(outdir_name))
    else:
        print('No outdir_name given, so will use default: {}'.format(AurnExtractor.DEFAULT_OUT_DIR))
        outdir_name = AurnExtractor.DEFAULT_OUT_DIR

    if args.statdir_name:
        statdir_name = args.statdir_name
        print('Using statdir_name: {}'.format(statdir_name))
    else:
        print('No statdir_name given, so will use default: {}'.format('./'))
        statdir_name = './'


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
        print('No verbose flag provided, so using default: {}'.format(str(AurnImputationTest.DEFAULT_VERBOSE)))
        verbose = AurnImputationTest.DEFAULT_VERBOSE

    extractor = AurnExtractor(metadata_filename=metadata_filename,
                               metadata_url=metadata_url,
                               out_dir=outdir_name, verbose=verbose)

    file_extracted = extractor._base_file_out.format(extractor.out_dir, '_'+outfile_suffix)

    processor = AurnImputationTest(metadata_filename=metadata_filename,
                                  metadata_url=metadata_url,
                                  out_dir=outdir_name, verbose=verbose,
                                  stat_dir=statdir_name)
    
    processor.impute_method_setup()
    
    processor.imputation_test(  in_file=file_extracted,
                        date_range=date_range,
                        outfile_suffix= outfile_suffix,
                        site_list=site_list,
                        species_list=species_list,
                        emep_filename=emep_filename,
                        min_years=min_years,
                        min_years_reference=min_years_ref,
                        data_lost=data_lost,
                        data_loss_position=data_loss_position,
                        save_to_csv=args.save_to_csv,
                        check_sites=args.check_sites)
