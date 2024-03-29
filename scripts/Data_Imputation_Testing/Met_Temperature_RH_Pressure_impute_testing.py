"""
Created on Tue Sep 22 09:16:40 2020

Script for testing the imputation of Met data.


@author:    Douglas Lowe and Ann Gledson
            Research IT, University of Manchester
"""

import argparse
import sys
sys.path.append("../..")
from pathlib import Path

#import readline # optional, will allow Up/Down/History in the console
#import code

from IPython import embed

from environmental_data_modules import MetImputationTest


if __name__ == '__main__':

    ### read arguments from the command line
    help_string = """Unified production script for MEDMI imputation testing."""

    parser = argparse.ArgumentParser(description=help_string)

    ## output directory/file names
    parser.add_argument("--out_dir", "-o", dest="out_dir", type=str,
                        help="output directory name. Default: {}".format(MetImputationTest.DEFAULT_OUT_DIR))
    parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
                        help="suffix to be appended to output file name. Default: '{}'".format(
                            MetImputationTest.DEFAULT_OUT_FILE_SUFFIX))
    parser.add_argument("--statdir_name", dest="statdir_name", type=str,
                        help="stat output directory name. Default: {}".format('./'))


    # input files
    parser.add_argument("--file_in", "-f", dest="file_in", type=str, help="name of input file. Required")
    parser.add_argument("--stations_filename", dest="stations_filename", type=str,
                        help="name of file containing list of Met stations. Required")

    ## Dates
    parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',
                        help="start and end dates. (array - first two values only). \
                            Expected date format: {} \n Default: [{}, {}]".format(
                            MetImputationTest.INPUT_DATE_FORMAT.replace('%', ''),
                            MetImputationTest.get_available_start(), MetImputationTest.get_available_end()))

    ## Calculation parameters
    parser.add_argument("--min_years", "-y", type=float, help="minimum number of years of data that a site must have")
    parser.add_argument("--min_years_ref", "-u", type=float, help="minimum number of years of data for any site that \
        we are going to use as a reference site later (if less than min_years, will set to min_years)")
    parser.add_argument("--ref_num_stations", "-r", type=int, help="number of stations to be used for imputation")
    parser.add_argument("--min_temp", "-m", type=int, help="Minimum temperature to be used (lower are ignored) \
                    (Default: {})".format(str(MetImputationTest.DEFAULT_MIN_TEMPERATURE)))
    parser.add_argument("--exclude_sites", "-x", metavar='X', dest="exclude_sites", type=str, nargs='+',
                        help="the measurement sites to be excluded. Default is zero sites.")

    parser.add_argument("--sites", metavar='S', dest="sites", type=str, nargs='+', help="the measurement sites \
        to be processed. Default is messed up at the moment, so don't use it!")

    # Imput values
    parser.add_argument("--impute_values", "-i", dest="impute_values", action='store_true',
                        help="impute missing values (default).")
    parser.add_argument("--no_impute_values", dest="impute_values", action='store_false',
                        help="don't impute missing values.")
    parser.set_defaults(impute_values=True)

    ## Printed outputs
    # Print statistics
    parser.add_argument("--print_stats", "-p", dest="print_stats", action='store_true',
                        help="print statistics (default).")
    parser.add_argument("--no_print_stats", dest="print_stats", action='store_false',
                        help="don't print statistics.")
    parser.set_defaults(impute_values=True)

    # Log verbose-ness
    parser.add_argument("--verbose", "-v", type=int,
                        help="Level of output for debugging (Default: {} (0 = no verbose output))".format(
                            str(MetImputationTest.DEFAULT_VERBOSE)))

    # Imputation test controls
    parser.add_argument("--data_lost", type=float, help="fraction of total period to remove for imputation tests")
    parser.add_argument("--data_loss_position", type=str, help="position to lose data: start, middle, end, random")
    parser.add_argument("--check_sites",dest="check_sites",action='store_true',help="list sites suitable for \
        tests then exit")
    parser.set_defaults(check_sites=False)

    # control interactive analysis work
    parser.add_argument("--interactive", dest="interactive", action='store_true', help="enter ipython interactive prompt for data analysis")
    parser.set_defaults(interactive=False)

    ### Process inputs
    args = parser.parse_args()

    if args.data_lost:
        data_lost = args.data_lost
        print('Removing {} fraction of data for imputation tests'.format(data_lost))
    else:
        data_lost = MetImputationTest.DEFAULT_DATA_LOST
        print('Default: Removing {} fraction of data for imputation tests'.format(data_lost))

    if args.data_loss_position:
        data_loss_position = args.data_loss_position
        print('Using method "{}" for selecting data to lose'.format(data_loss_position))
    else:
        data_loss_position = MetImputationTest.DEFAULT_DATA_LOSS_POSITION
        print('Default: Using method "{}" for selecting data to lose'.format(data_loss_position))

    if args.file_in:
        file_in = Path(args.file_in)
    else:
        raise ValueError('Error: No input filename provided')

    if args.stations_filename:
        stations_filename = Path(args.stations_filename)
    else:
        raise ValueError('Error: No stations filename provided')

    if args.out_dir:
        out_dir = args.out_dir
        print('Using outdir_name: {}'.format(out_dir))
    else:
        print('No out_dir given, so will use default: {}'.format(MetImputationTest.DEFAULT_OUT_DIR))
        outdir_name = MetImputationTest.DEFAULT_OUT_DIR

    if args.outfile_suffix:
        outfile_suffix = args.outfile_suffix
        print('Using outfile_suffix: {}'.format(outfile_suffix))
    else:
        print('No outfile_suffix provided, so no suffix added')
        outfile_suffix = ''

    if args.statdir_name:
        statdir_name = args.statdir_name
        print('Using statdir_name: {}'.format(statdir_name))
    else:
        print('No statdir_name given, so will use default: {}'.format('./'))
        statdir_name = './'

    if args.date_range:
        if len(args.date_range) >= 2:
            if len(args.date_range) > 2:
                print('Warning: You have input more than 2 dates, only the first 2 will be used.')
            date_range = args.date_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 dates from input --date_range: {}'.format(str(args.date_range)))
        print('Using date range: [{}]'.format(','.join(date_range)))
    else:
        print('No date_range provided, so using default: [{}]'.format(','.join(MetImputationTest.get_available_dates())))
        date_range = MetImputationTest.get_available_dates()

    if args.min_years:
        min_years = args.min_years
        print('Min years (minimum number of years of data that a site must have):', min_years)
    else:
        print('No min_years provided, so using default: '.format(MetImputationTest.DEFAULT_MIN_YEARS))
        min_years = MetImputationTest.DEFAULT_MIN_YEARS

    if args.min_years_ref:
        min_years_ref = max(args.min_years_ref, min_years)
        print('Reference number of years (minimum number of years of data for any site that we are going to use as a \
            reference site later; this cannot be less than min_years):', min_years_ref)
    else:
        print('No min_years_ref provided, so using default: {}'.format(MetImputationTest.DEFAULT_MIN_YEARS_REFERENCE))
        min_years_ref = MetImputationTest.DEFAULT_MIN_YEARS_REFERENCE

    if args.ref_num_stations:
        reference_num_stations = args.ref_num_stations
        print('Reference number of stations :', reference_num_stations)
    else:
        print('No useful_num_stations provided, so using default: {}'.format(
            MetImputationTest.DEFAULT_REFERENCE_NUMBER_STATIONS))
        useful_num_stations = MetImputationTest.DEFAULT_REFERENCE_NUMBER_STATIONS

    if args.min_temp:
        min_temperature = args.min_temp
        print('Min temperature (minimum tempoerature reading (or will be ignored):', min_temperature)
    else:
        print('No min_temperature provided, so using default: '.format(MetImputationTest.DEFAULT_MIN_TEMPERATURE))
        min_temperature = MetImputationTest.DEFAULT_MIN_TEMPERATURE

    if args.sites:
        site_list = args.sites
        print('Site list: {}'.format(site_list))
    else:
        print('No sites provided, so using all available sites in metadata file')
        site_list = None

    if args.exclude_sites:
        exclude_site_list = args.exclude_sites
        print('Excluding sites: [{}]'.format(','.join(exclude_site_list)))
    else:
        print('No exclude_sites provided, so using all')
        exclude_site_list = []

    print('Impute values: {}'.format(args.impute_values))
    print('Print stats: {}'.format(args.print_stats))

    if args.verbose:
        verbose = max(args.verbose, 0)
        print('output verbose level: {}'.format(verbose))
    else:
        print('No verbose flag provided, so using default: {}'.format(str(MetImputationTest.DEFAULT_VERBOSE)))
        verbose = MetImputationTest.DEFAULT_VERBOSE

    post_processor = MetImputationTest(out_dir, station_data_filename=stations_filename, 
                     stat_dir=statdir_name, verbose=verbose)
    post_processor.impute_method_setup()

    post_processor.imputation_test(file_in, outfile_suffix=outfile_suffix, date_range=date_range,
                                exclude_site_list=exclude_site_list,
                                min_temperature=min_temperature, reference_num_stations=reference_num_stations,
                                min_years=min_years, min_years_reference=min_years_ref,
                                print_stats=args.print_stats,
                                data_lost=data_lost,data_loss_position=data_loss_position,
                                check_sites=args.check_sites,site_list=site_list)
