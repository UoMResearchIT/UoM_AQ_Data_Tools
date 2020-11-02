import argparse
import sys
sys.path.append("..")
from pathlib import Path

from environmental_data_modules import MetPostProcessor


if __name__ == '__main__':

    ### read arguments from the command line
    help_string = "***Unified production script for MEDMI data processing. ***\
        Loads and cleans the dataset by:\
        a) find duplicated readings \
            i)  filtering out the METAR data by lack of pressure reading\
            ii) retaining the first value where there's no difference in presence of pressure reading\
        b) removing stations identified as unwanted\
            i) station 117 is on top of a mountain in the Cairngorms - RH readings are suspect,\
                    and as it is unlikely to be useful comparison with participant data, we will remove it\
        c) find and remove the synoptic spot readings"

    parser = argparse.ArgumentParser(description=help_string)

    ## output directory/file names
    parser.add_argument("--out_dir", "-o", dest="out_dir", type=str,
                        help="output directory name. Default: {}".format(MetPostProcessor.DEFAULT_OUT_DIR))
    parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
                        help="suffix to be appended to output file name. Default: {}".format(
                            MetPostProcessor.DEFAULT_OUT_FILE_SUFFIX))

    # input files
    parser.add_argument("--file_in", "-f", dest="file_in", type=str, help="name of input file. Required")
    parser.add_argument("--stations_filename", dest="stations_filename", type=str,
                        help="name of file containing list of Met stations. Required")

    ## Dates
    parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',
                        help="start and end dates. Expected date format: {}. \
                            (array - first two values only). Default: [{}]".format(
                            MetPostProcessor.DEFAULT_DATE_RANGE_FORMAT,
                            ','.join(MetPostProcessor.DEFAULT_DATE_RANGE)))

    ## Calculation parameters
    parser.add_argument("--min_years", "-y", type=float, help="minimum number of years of data that a site must have")
    parser.add_argument("--ref_num_years", "-u", type=float, help="minimum number of years of data for any site that \
        we are going to use as a reference site later (if less than min_years, will set to min_years)")
    parser.add_argument("--ref_num_stations", "-r", type=int, help="number of stations to be used for imputation")
    parser.add_argument("--min_temp", "-m", type=int, help="Minimum temperature to be used (lower are ignored) \
                    (Default: {})".format(str(MetPostProcessor.DEFAULT_MIN_TEMPERATURE)))
    parser.add_argument("--exclude_sites", "-x", metavar='X', dest="exclude_sites", type=str, nargs='+',
                        help="the measurement sites to be excluded. Default is zero sites.")

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
                            str(MetPostProcessor.DEFAULT_VERBOSE)))

    ### Process inputs
    args = parser.parse_args()

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
        print('No out_dir given, so will use default: {}'.format(MetPostProcessor.DEFAULT_OUT_DIR))
        outdir_name = MetPostProcessor.DEFAULT_OUT_DIR

    if args.outfile_suffix:
        outfile_suffix = args.outfile_suffix
        print('Using outfile_suffix: {}'.format(outfile_suffix))
    else:
        print('No outfile_suffix provided, so no suffix added')
        outfile_suffix = ''

    if args.date_range:
        if len(args.date_range) >= 2:
            if len(args.date_range) > 2:
                print('Warning: You have input more than 2 dates, only the first 2 will be used.')
            date_range = args.date_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 dates from input --date_range: {}'.format(str(args.date_range)))
        print('Using date range: [{}]'.format(','.join(date_range)))
    else:
        print('No date_range provided, so using default: [{}]'.format(','.join(MetPostProcessor.DEFAULT_DATE_RANGE)))
        date_range = MetPostProcessor.DEFAULT_DATE_RANGE

    if args.min_years:
        min_years = args.min_years
        print('Min years (minimum number of years of data that a site must have):', min_years)
    else:
        print('No min_years provided, so using default: '.format(MetPostProcessor.DEFAULT_MIN_YEARS))
        min_years = MetPostProcessor.DEFAULT_MIN_YEARS

    if args.ref_num_years:
        reference_num_years = max(args.ref_num_years, min_years)
        print('Reference number of years (minimum number of years of data for any site that we are going to use as a \
            reference site later; this cannot be less than min_years):', reference_num_years)
    else:
        print('No useful_num_years provided, so using default: {}'.format(MetPostProcessor.DEFAULT_REFERENCE_NUM_YEARS))
        useful_num_years = MetPostProcessor.DEFAULT_REFERENCE_NUM_YEARS

    if args.ref_num_stations:
        reference_num_stations = args.ref_num_stations
        print('Reference number of stations :', reference_num_stations)
    else:
        print('No useful_num_stations provided, so using default: {}'.format(
            MetPostProcessor.DEFAULT_REFERENCE_NUMBER_STATIONS))
        useful_num_stations = MetPostProcessor.DEFAULT_REFERENCE_NUMBER_STATIONS

    if args.min_temp:
        min_temperature = args.min_temp
        print('Min temperature (minimum tempoerature reading (or will be ignored):', min_temperature)
    else:
        print('No min_temperature provided, so using default: '.format(MetPostProcessor.DEFAULT_MIN_TEMPERATURE))
        min_temperature = MetPostProcessor.DEFAULT_MIN_TEMPERATURE

    if args.exclude_sites:
        exclude_site_list = args.exclude_sites
        print('Excluding sites: [{}]'.format(','.join(exclude_site_list)))
    else:
        print('No exclude_sites provided, so using all')
        exclude_site_list = []

    print('Impute values: {}'.format(args.impute_values))
    print('Print stats: {}'.format(args.impute_values))

    if args.verbose:
        verbose = max(args.verbose, 0)
        print('output verbose level: {}'.format(verbose))
    else:
        print('No verbose flag provided, so using default: {}'.format(str(MetPostProcessor.DEFAULT_VERBOSE)))
        verbose = MetPostProcessor.DEFAULT_VERBOSE

    post_processor = MetPostProcessor(out_dir, stations_filename=stations_filename, verbose=verbose)
    post_processor.post_process(file_in, outfile_suffix=outfile_suffix, date_range=date_range,
                                exclude_site_list=exclude_site_list,
                                min_temperature=min_temperature, reference_num_stations=reference_num_stations,
                                min_years=min_years, reference_num_years=reference_num_years,
                                impute_data=args.impute_values, print_stats=args.print_stats)