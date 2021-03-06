"""
Created on Tue Sep 22 09:16:40 2020

Unified production script for MEDMI data extraction. Must be run on the MEDMI server.  To use these script(s),
it is neccessary to log onto the MEDMI ssh server. Details of how to do this can be obtained from the
MEDMI website: https://www.data-mashup.org.uk/contact-us/
and following the instructions for: connect to the server using an SSH client.
Alternatively, email: health at metoffice dot gov dot uk

Copy this script to the ssh server. Then run it using python (simply python [script name] [params/options])
and use the --help command to get the full set of parameters/options.

    This will:
    1) Extract weathr and pollen data from the MEDMI server, given the inputs:
        a) Location: expressed as latitude and longitude ranges
        b) Dates: expressed as a date range (start to end)
        c) Required measurements
    ToDo: Doug ++?

@author:    Douglas Lowe and Ann Gledson
            Research IT, University of Manchester
"""

import argparse
import sys
sys.path.append("../..")
from environmental_data_modules import MetExtractor

if __name__ == '__main__':

    ### general help text
    parser = argparse.ArgumentParser(
        description="*** A script for automated downloading of MET data for a given region and date range. ***")

    AVAILABLE_MEASUREMENTS = MetExtractor.get_available_measurements()

    ### parameters

    # Measurements
    parser.add_argument("--measurements", "-m", metavar='M', type=str, nargs='+', help="measurements to be extracted. \
            Must be in (and defaults to) {}".format('[' + ", ".join([measure for measure in AVAILABLE_MEASUREMENTS]) + ']'))

    parser.add_argument("--extra_measurements", "-x", dest='extra_measurements', action='store_true', \
                        help="request extra measurements to be  extracted along with each of the main measurements.")

    # output directory/file names
    parser.add_argument("--outdir_name", "-o", dest="outdir_name", type=str,
                        help="output directory name. Default: {}".format(MetExtractor.DEFAULT_OUT_DIR))
    parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
                        help="suffix to be appended to output file name. Default: {}".format(MetExtractor.DEFAULT_OUT_FILE_SUFFIX))

    # Dates
    parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',
                        help="start and end dates. (array - first two values only). \
                            Expected date format: {} \n Default: [{}, {}]".format(
                            MetExtractor.INPUT_DATE_FORMAT.replace('%', ''),
                            MetExtractor.get_available_start(), MetExtractor.get_available_end()))


    # Latitude / longitude
    parser.add_argument("--latitude_range", "-t", dest="latitude_range", type=float, nargs='+',
                        help="start and end latitude range (array - first two values only). \
                            Default: {}".format(str(MetExtractor.UK_LATITUDES)[1:-1].replace(',','')))
    parser.add_argument("--longitude_range", "-n", dest="longitude_range", type=float, nargs='+',
                        help="start and end longitude range (array - first two values only). \
                            Default: {}".format(str(MetExtractor.UK_LONGITUDES)[1:-1].replace(',','')))

    # Log verbose-ness
    parser.add_argument("--verbose", "-v", type=int,
                        help="Level of output for debugging (Default: {} (0 = verbose output))".format(str(MetExtractor.DEFAULT_VERBOSE)))

    ### Process Inputs
    args = parser.parse_args()

    if args.measurements:
        for measurement in args.measurements:
            if measurement not in AVAILABLE_MEASUREMENTS:
                raise ValueError('Unknown measurement: {}. Allowed measurements: {}'.format(measurement, AVAILABLE_MEASUREMENTS))
        # Remove duplicate measurement names.
        measurements = set(args.measurements)
        print('Using measurements: {}'.format(measurements))
    else:
        print('No measurements provided, so using default: ', '[' + ", ".join([measure for measure in AVAILABLE_MEASUREMENTS]) + ']')
        measurements = AVAILABLE_MEASUREMENTS

    if args.extra_measurements is True:
        extra_measurements = True
        print('Including extra measurements')
    else:
        print('Not including extra measurements')
        extra_measurements = False

    if args.outdir_name:
        outdir_name = args.outdir_name
        print('Using outdir_name: {}'.format(outdir_name))
    else:
        print('No outdir_name given, so will use default: {}'.format(MetExtractor.DEFAULT_OUT_DIR))
        outdir_name = MetExtractor.DEFAULT_OUT_DIR

    if args.outfile_suffix:
        outfile_suffix = args.outfile_suffix
        print('Using outfile_suffix: {}'.format(outfile_suffix))
    else:
        print('No outfile_suffix provided, so using default: {}'.format(str(MetExtractor.DEFAULT_OUT_FILE_SUFFIX)))
        outfile_suffix = MetExtractor.DEFAULT_OUT_FILE_SUFFIX

    if args.date_range:
        #'2017-01-01_00 2019-06-30_23'
        if len(args.date_range) >= 2:
            if len(args.date_range) > 2:
                print('Warning: You have input more than 2 dates, only the first 2 will be used.')
            date_range = args.date_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 dates from input --date_range: {}'.format(str(args.date_range)))
        print('Using date range: [{}]'.format(','.join(date_range)))
    else:
        print('No date_range provided, so using default: [{}, {}]'.format(MetExtractor.get_available_start(),
                                                                         MetExtractor.get_available_end()))
        date_range = [MetExtractor.get_available_start(), MetExtractor.get_available_end()]

    if args.latitude_range:
        if len(args.latitude_range) >= 2:
            if len(args.latitude_range) > 2:
                print('Warning: You have input more than 2 latitudes, only the first 2 will be used.')
            latitude_range = args.latitude_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 values from input --latitude_range: {}'.format(str(args.latitude_range)))
        print('Using latitude range: {}'.format(str(latitude_range)[1:-1].replace(',','')))
    else:
        print('No latitude_range provided, so using default: {}'.format(str(MetExtractor.UK_LATITUDES)[1:-1].replace(',','')))
        latitude_range = MetExtractor.UK_LATITUDES

    if args.longitude_range:
        if len(args.longitude_range) >= 2:
            if len(args.longitude_range) > 2:
                print('Warning: You have input more than 2 longitudes, only the first 2 will be used.')
            longitude_range = args.longitude_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 values from input --longitude_range: {}'.format(str(args.longitude_range)))
        print('Using longitude range: {}'.format(str(longitude_range)[1:-1].replace(',','')))
    else:
        print('No longitude_range provided, so using default: {}'.format(str(MetExtractor.UK_LONGITUDES)[1:-1].replace(',','')))
        longitude_range = MetExtractor.UK_LONGITUDES

    if args.verbose:
        verbose = max(args.verbose, 0)
        print('output verbose level: {}'.format(verbose))
    else:
        print('No verbose flag provided, so using default: {}'.format(str(MetExtractor.DEFAULT_VERBOSE)))
        verbose = MetExtractor.DEFAULT_VERBOSE

    ### Perform data extraction

    # Perform extractions
    for measurement in measurements:
        class_ = MetExtractor.get_class_from_measurement_name(measurement)
        met_extractor = class_(outdir_name, verbose)
        met_extractor.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)