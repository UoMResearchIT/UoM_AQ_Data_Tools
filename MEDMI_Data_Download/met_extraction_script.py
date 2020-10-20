#from abc import ABCMeta, abstractmethod
from medmi_database import Dataset
from cmath import polar
import argparse

import os, errno


class MetExtractor(object):
    #__metaclass__ = ABCMeta

    VERBOSE = 0
    DEFAULT_OUT_DIR_PREFIX = 'data_'
    DEFAULT_OUT_FILE_SUFFIX = '2016-2019'
    DEFAULT_DATE_RANGE = ['2016-1-1 0', '2019-12-31 23']
    DEFAULT_LATITUDES = [48, 60]
    DEFAULT_LONGITUDES = [-11, 3]
    DEFAULT_VERBOSE = 0
    ADDITIONAL_VALS_KEY = 'Additional field values'


    def __init__(self, dir_name=DEFAULT_OUT_DIR_PREFIX):
        self.out_dir = dir_name
        self.cols_base = 'date,siteID,{}\n'
        self.outfile_suffix = MetExtractor.DEFAULT_OUT_FILE_SUFFIX

    def _create_dict_base(self, date_range, latitude_range, longitude_range):
        return {'Time range': date_range, 'Latitude range': latitude_range, 'Longitude range': longitude_range}

    def _extract_data(self, dict, settings, extra_datasets=[]):
        datadata = Dataset(dict)
        datadata.default()
        for ds in extra_datasets:
            datadata.add(ds)

        with open(settings['fname'], 'w') as dfile:
            dfile.write(settings['headstring'])
            dfile.write(settings['columnstring'])

            for data in datadata.values():
                d_date = data['Time']
                d_siteid = data['Site identifier']
                d_val = data['Value']

                if MetExtractor.ADDITIONAL_VALS_KEY in data:
                    d_extra = data[MetExtractor.ADDITIONAL_VALS_KEY]
                else:
                    d_extra = []

                dfile.write('{}, {}, {}'.format(d_date, d_siteid, d_val))
                for d_ev in d_extra:
                    dfile.write(', {}'.format(d_ev))
                dfile.write('\n')


class MetExtractorRain(MetExtractor):
    aliases_lowercase = ['rain', 'rainfall', 'precipitation']

    def __init__(self, outfile_dir):
        super(MetExtractorRain, self).__init__(outfile_dir)
        self.source_reference = 'midas.rain_drnl_ob.prcp_amt 1'
        self.head_string = 'Rain gauge daily data for date range: {} to {}\n'

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        column_string = self.cols_base.format('rain')
        filename = '{}/rain_{}.csv'.format(self.out_dir, outfile_suffix)
        print('extracting {}'.format(self.head_string.format(date_range[0], date_range[1])))
        dict = super(MetExtractorRain, self)._create_dict_base(date_range, latitude_range, longitude_range)
        dict.update({'Source reference': self.source_reference})
        settings = {'fname': filename,
                    'headstring': self.head_string.format(date_range[0], date_range[1]),
                    'columnstring': column_string}
        self._extract_data(dict, settings)


class MetExtractorTemp(MetExtractor):
    aliases_lowercase = ['temp', 'temperature']

    def __init__(self, outfile_dir):
        super(MetExtractorTemp, self).__init__(outfile_dir)
        self.source_reference = 'midas.weather_hrly_ob.air_temperature'
        self.head_string = 'Temperature, Relative Humidity, Station Pressure, and Wet Bulb Temperature hourly \
                                    data for date range: {} to {}\n'

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        column_string = self.cols_base.format('temperature,rh,pressure,wbtemp')
        filename = '{}/temp_rh_press_wbulb_{}.csv'.format(new_dir, outfile_suffix)
        extra_datasets = ['midas.weather_hrly_ob.rltv_hum', 'midas.weather_hrly_ob.stn_pres',
                          'midas.weather_hrly_ob.dewpoint']
        print('extracting {}'.format(self.head_string.format(date_range[0], date_range[1])))
        dict = self._create_dict_base(date_range, latitude_range, longitude_range)
        dict.update({'Source reference': self.source_reference})
        settings = {'fname': filename,
                    'headstring': self.head_string.format(date_range[0], date_range[1]),
                    'columnstring': column_string}
        self._extract_data(dict, settings, extra_datasets)



class MetExtractorWind(MetExtractor):

    def __init__(self, outfile_dir):
        super(MetExtractorWind, self).__init__(outfile_dir)
        self.source_reference = 'midas.wind_mean_ob.mean_wind_speed 1'
        self.complex_wind_type = True
        self.head_string = 'Wind speed and direction hourly data for date range: {} to {}\n'

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        self.outfile_suffix = outfile_suffix
        self.filename = '{}/wind_{}.csv'.format(new_dir, outfile_suffix)
        print('extracting {}'.format(self.head_string.format(date_range[0], date_range[1])))
        dict = super(MetExtractorWind, self)._create_dict_base(date_range, latitude_range, longitude_range)
        dict.update({'Source reference': self.source_reference, 'Complex wind type': self.complex_wind_type})

        settings = {    'fname': self.filename,
                        'headstring': self.head_string.format(date_range[0], date_range[1]),
                        'columnstring': self.cols_base.format('windspeed,winddir')}

        self._extract_data(dict, settings)

    # Overrides base class method
    def _extract_data(self, dict, settings):
        datadata = Dataset(dict)
        datadata.default()

        with open(settings['fname'], 'w') as dfile:
            dfile.write(settings['headstring'])
            dfile.write(settings['columnstring'])
            for data in datadata.values():
                d_date = data['Time']
                d_siteid = data['Site identifier']
                (d_wspd, d_wdir) = polar(data['Value'])
                d_wdir = d_wdir * 57.29577951308232
                if d_wdir < 0: d_wdir = d_wdir + 360.
                dfile.write('{}, {}, {}, {}\n'.format(d_date, d_siteid, d_wspd, d_wdir))


class MetExtractorPollen(MetExtractor):

    def __init__(self, outfile_dir):
        super(MetExtractorPollen, self).__init__(outfile_dir)
        self.source_reference_base = 'midas.pollen_drnl_ob.{}'
        self.pollen_species = ['alnus', 'ambrosia', 'artemisia', 'betula', 'corylus', 'fraxinus', 'platanus', 'poaceae',
                          'quercus', 'salix', 'ulmus', 'urtica']
        self.head_string = '{} pollen daily count for date range: {} to {}\n'

    def extract_all_datasets(self, date_range, latitude_range, longitude_range, outfile_suffix):
        self.filename = '{}/pollen_{}_{}.csv'.format(new_dir, outfile_suffix, '{}')
        dict_base = super(MetExtractorPollen, self)._create_dict_base(date_range, latitude_range, longitude_range)
        for pspc in self.pollen_species:
            print('working on pollen species: {}'.format(pspc))
            fname = self.filename.format(pspc)
            dict = dict_base.copy()
            dict['Source reference'] = self.source_reference_base.format(pspc)
            settings = {'fname': fname,
                        'headstring': self.head_string.format(pspc, dict_base['Time range'][0], dict_base['Time range'][1]),
                        'columnstring': self.cols_base.format(pspc)}
            super(MetExtractorPollen, self)._extract_data(dict, settings)



def create_directory(dir_name):
    try:
        os.makedirs(dir_name)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


if __name__ == '__main__':

    ### general help text
    parser = argparse.ArgumentParser(
        description="*** A script for automated downloading of MET data for a given date range. ***")

    AVAILABLE_MEASUREMENTS = ['pollen', 'rain', 'temp',  'wind']

    ### parameters

    # Measurement
    parser.add_argument("--measurements", "-m", metavar='M', type=str, nargs='+', help="measurements to be extracted. \
            Must be in (and defaults to) {}".format('[' + ", ".join([measure for measure in AVAILABLE_MEASUREMENTS]) + ']'))
    #parser.add_argument("--measurements", "-m", dest="measurement", type=str, nargs='+',
    #                    help="measurements to be extracted. Default: all")

    # output directory/file names
    parser.add_argument("--outdir_prefix", "-o", dest="outdir_prefix", type=str,
                        help="prefix to be appended on output directory name. Default: {}".format(MetExtractor.DEFAULT_OUT_DIR_PREFIX))
    parser.add_argument("--outfile_suffix", "-s", dest="outfile_suffix", type=str,
                        help="suffix to be appended to output file name. Default: {}".format(MetExtractor.DEFAULT_OUT_FILE_SUFFIX))

    # Dates
    parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',  help="start and end dates \
                        (array - first two values only). Default: {}".format(str(MetExtractor.DEFAULT_DATE_RANGE)[1:-1].replace(',','')))


    # Latitude / longitude
    parser.add_argument("--latitude_range", "-t", dest="latitude_range", type=int, nargs='+',
                        help="start and end latitude range (array - first two values only). \
                            Default: {}".format(str(MetExtractor.DEFAULT_LATITUDES)[1:-1].replace(',','')))
    parser.add_argument("--longitude_range", "-n", dest="longitude_range", type=int, nargs='+',
                        help="start and end longitude range (array - first two values only). \
                            Default: {}".format(str(MetExtractor.DEFAULT_LONGITUDES)[1:-1].replace(',','')))

    # Log verbose-ness
    parser.add_argument("--verbose", "-v", type=int,
                        help="Level of output for debugging (Default: {} (0 = verbose output))".format(str(MetExtractor.DEFAULT_VERBOSE)))

    ### Process Inputs
    args = parser.parse_args()

    if args.measurements:
        measurements = args.measurements
        print('Using measurements: {}'.format(measurements))
    else:
        print('No measurements provided, so using default: ', '[' + ", ".join([measure for measure in AVAILABLE_MEASUREMENTS]) + ']')
        measurements = AVAILABLE_MEASUREMENTS

    if args.outdir_prefix:
        outdir_prefix = args.outdir_prefix
        print('Using outdir_prefix: {}'.format(outdir_prefix))
    else:
        print('No outdir_prefix given, so will use default: {}'.format(MetExtractor.DEFAULT_OUT_DIR_PREFIX))
        outdir_prefix = MetExtractor.DEFAULT_OUT_DIR_PREFIX

    if args.outfile_suffix:
        outfile_suffix = args.outfile_suffix
        print('Using outfile_suffix: {}'.format(outfile_suffix))
    else:
        print('No outfile_suffix provided, so using default: {}'.format(str(MetExtractor.DEFAULT_OUT_FILE_SUFFIX)))
        outfile_suffix = MetExtractor.DEFAULT_OUT_FILE_SUFFIX

    if args.date_range:
        if len(args.date_range) >= 2:
            date_range = args.date_range[0:2]
        else:
            raise ValueError(
                'Unable to obtain 2 dates from input --date_range: {}'.format(str(args.date_range)))
        print('Using date range: {}'.format(str(date_range)[1:-1].replace(',','')))
    else:
        print('No date_range provided, so using default: {}'.format(str(MetExtractor.DEFAULT_DATE_RANGE)[1:-1].replace(',','')))
        date_range = MetExtractor.DEFAULT_DATE_RANGE

    if args.latitude_range:
        if len(args.latitude_range) >= 2:
            latitude_range = args.latitude_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 values from input --latitude_range: {}'.format(str(args.latitude_range)))
        print('Using latitude range: {}'.format(str(latitude_range)[1:-1].replace(',','')))
    else:
        print('No latitude_range provided, so using default: {}'.format(str(MetExtractor.DEFAULT_LATITUDES)[1:-1].replace(',','')))
        latitude_range = MetExtractor.DEFAULT_LATITUDES

    if args.longitude_range:
        if len(args.longitude_range) >= 2:
            longitude_range = args.longitude_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 values from input --longitude_range: {}'.format(str(args.longitude_range)))
        print('Using longitude range: {}'.format(str(longitude_range)[1:-1].replace(',','')))
    else:
        print('No longitude_range provided, so using default: {}'.format(str(MetExtractor.DEFAULT_LONGITUDES)[1:-1].replace(',','')))
        longitude_range = MetExtractor.DEFAULT_LONGITUDES

    if args.verbose:
        VERBOSE = max(args.verbose, 0)
        print('verbose: ', VERBOSE)
    else:
        print('No verbose flag provided, so using default: {}'.format(str(MetExtractor.DEFAULT_VERBOSE)))
        VERBOSE = MetExtractor.DEFAULT_VERBOSE


    # Check directory name is valid and create directory
    new_dir = '{}met'.format(outdir_prefix)
    try:
        print('Creating directory: {}, unless it already exists.'.format(new_dir))
        create_directory(new_dir)
    except:
        raise ValueError('Unable to create directory/file: {}'.format(new_dir))


    ### Prepare inputs and perform data extraction

    # Create necessary extractors
    if 'rain' in measurements:
        met_extractor_rain = MetExtractorRain(new_dir)
        met_extractor_rain.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
    if 'temp' in measurements:
        met_extractor_temp = MetExtractorTemp(new_dir)
        met_extractor_temp.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
    if 'wind' in measurements:
        met_extractor_wind = MetExtractorWind(new_dir)
        met_extractor_wind.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
    if 'pollen' in measurements:
        met_extractor_pollen = MetExtractorPollen(new_dir)
        met_extractor_pollen.extract_all_datasets(date_range, latitude_range, longitude_range, outfile_suffix)



