from abc import ABCMeta, abstractmethod
from medmi_database import Dataset
from cmath import polar
import argparse
import json

import os, errno


class MetExtractor:
    __metaclass__ = ABCMeta
    DEFAULT_OUT_FILE_SUFFIX = ''
    DEFAULT_OUT_DIR = 'met_extracted_data'
    DEFAULT_DATE_RANGE = ['2016-1-1 0', '2019-12-31 23']
    UK_LATITUDES = [48, 60]
    UK_LONGITUDES = [-11, 3]
    DEFAULT_VERBOSE = 0
    DEFAULT_COLS_BASE = ['date','siteID']
    ADDITIONAL_VALS_KEY = 'Additional field values'
    DEFAULT_ADD_EXTRA_MEASUREMENTS = False


    def __init__(self, dir_name=DEFAULT_OUT_DIR, verbose=DEFAULT_VERBOSE):
        self._out_dir = dir_name
        self._verbose = verbose
        self._cols_base = MetExtractor.DEFAULT_COLS_BASE
        self._latitude_range = MetExtractor.UK_LATITUDES
        self._longitude_range = MetExtractor.UK_LONGITUDES
        self._date_range = MetExtractor.DEFAULT_DATE_RANGE
        self._filename = None
        self._headstring = None
        self._extra_datasets = []


    @staticmethod
    def get_source_ref_from_name(name):
        for sub_class in MetExtractor.all_subclasses(MetExtractor):
            if hasattr(sub_class, 'MEASUREMENT_NAME') and sub_class.MEASUREMENT_NAME == name:
                return sub_class.SOURCE_REFERENCE

    @staticmethod
    def get_top_level_measurements():
        return [s.MEASUREMENT_NAME for s in MetExtractor.__subclasses__()]

    @staticmethod
    def all_subclasses(cls):
        return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in MetExtractor.all_subclasses(c)])


    def add_outfile_suffix(self, outfile_suffix):
        return self._filename.format(outfile_suffix)

    def set_date_range(self, date_range):
        self._date_range = date_range

    def get_date_range(self):
        return self._date_range

    def get_out_dir(self):
        return self._out_dir

    def set_region(self, latitude_range, longitude_range):
        self._latitude_range = latitude_range
        self._longitude_range = longitude_range


    def _get_extraction_dict(self):
        result = {
            'Time range': self._date_range,
            'Source reference': self.SOURCE_REFERENCE,
            'Latitude range': self._latitude_range,
            'Longitude range': self._longitude_range
        }

        return result

    def _get_settings(self, outfile_suffix=DEFAULT_OUT_FILE_SUFFIX):
        str_extra_datasets = str(self._extra_datasets).replace('\'','').replace(' ', '').replace(',', '-')[1:-1]
        if len(self._extra_datasets) > 0:
            filename = self._filename.format('_extras-{}'.format(str_extra_datasets), outfile_suffix)
            headstring = self._head_string.format(' and extra datasets: {}'.format(str(str_extra_datasets)),
                                                   date_range[0], date_range[1])
        else:
            filename = self._filename.format('', outfile_suffix)
            headstring = self._head_string.format('', date_range[0], date_range[1])

        result = {
            'fname': filename,
            'headstring': headstring,
            'columnstring': str(self._cols_base + self._cols_specific + self._extra_datasets)[1:-1] + '\n'
        }
        return result



    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets=False):
        self.set_region(latitude_range, longitude_range)
        self.set_date_range(date_range)

        if extract_extra_datasets:
            self._extra_datasets = self.ALLOWED_EXTRA_DATASETS
        else:
            self._extra_datasets = []\

        extraction_dict = self._get_extraction_dict()
        settings = self._get_settings(outfile_suffix)
        print('extracting {}'.format(settings['headstring'], date_range[0], date_range[1]))
        return self._extract_data(extraction_dict, settings)


    def _extract_data(self, extraction_dict, settings, save_to_file=True):
        if self._verbose >= 1:
            print('extracting data for {}'.format(self.MEASUREMENT_NAME))
        if self._verbose > 1:
            print('using extraction dict: {}'.format(json.dumps(extraction_dict)))
            print('using settings: {}'.format( json.dumps(settings)))
        datadata = self._perform_extraction(extraction_dict)
        if save_to_file:
            self._save_to_file(datadata, settings)
        return datadata

    def _perform_extraction(self, dict):
        datadata = Dataset(dict)
        datadata.default()

        for ds in self._extra_datasets:
            source_reference = self.__class__.get_source_ref_from_name(ds)
            if source_reference is not None:
                print('extracting extra-data: {}, {}'.format(ds, source_reference))
                datadata.add(source_reference)
            else:
                print('could not find individual source reference for {}'.format(ds))

        return datadata

    def _save_to_file(self, data_result, settings):
        print('saving to file: {}'.format(settings['fname']))
        with open(settings['fname'], 'w') as dfile:
            dfile.write(settings['headstring'])
            dfile.write(settings['columnstring'])

            for data in data_result.values():
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
        if self._verbose > 0:
            print('saved')





class MetExtractorRain(MetExtractor):
    SOURCE_REFERENCE = 'midas.rain_drnl_ob.prcp_amt 1'
    MEASUREMENT_NAME = 'rain'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorRain, self).__init__(out_dir, verbose)
        self._head_string = 'Rain gauge daily data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorRain.MEASUREMENT_NAME]
        self._filename = '{}/rain{}{}.csv'.format(self._out_dir, '{}', '{}')



class MetExtractorTemp(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.air_temperature'
    MEASUREMENT_NAME = 'temp'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'pressure', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorTemp, self).__init__(out_dir, verbose)
        self._head_string = 'Temperature data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorTemp.MEASUREMENT_NAME]
        self._filename = '{}/temp{}{}.csv'.format(self._out_dir, '{}', '{}')


class MetExtractorRelativeHumidity(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.rltv_hum'
    MEASUREMENT_NAME = 'rel_hum'
    ALLOWED_EXTRA_DATASETS = ['temp', 'pressure', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorRelativeHumidity, self).__init__(out_dir, verbose)
        self._head_string = 'Relative Humidity data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorRelativeHumidity.MEASUREMENT_NAME]
        self._filename = '{}/rel_hum{}{}.csv'.format(self._out_dir, '{}', '{}')


class MetExtractorStationPressure(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.stn_pres'
    MEASUREMENT_NAME = 'pressure'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'temp', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorStationPressure, self).__init__(out_dir, verbose)
        self._head_string = 'Station Pressure data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorStationPressure.MEASUREMENT_NAME]
        self._filename = '{}/pressure{}{}.csv'.format(self._out_dir, '{}', '{}')



class MetExtractorDewpoint(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.dewpoint'
    MEASUREMENT_NAME = 'dewpoint'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'pressure', 'temp']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorDewpoint, self).__init__(out_dir, verbose)
        self._head_string = 'Dewpoint hourly data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorDewpoint.MEASUREMENT_NAME]
        self._filename = '{}/dewpoint{}{}.csv'.format(self._out_dir, '{}', '{}')


class MetExtractorWind(MetExtractor):
    SOURCE_REFERENCE = 'midas.wind_mean_ob.mean_wind_speed 1'
    MEASUREMENT_NAME = 'wind'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorWind, self).__init__(out_dir, verbose)
        self._complex_wind_type = True
        self._head_string = 'Wind speed/direction hourly data{} for date range: {} to {}\n'
        self._cols_specific = ['windspeed','winddir']
        self._filename = '{}/wind{}{}.csv'.format(out_dir, '{}', '{}')

    def _get_extraction_dict(self):
        result = super(MetExtractorWind, self)._get_extraction_dict()
        result['Complex wind type'] = self._complex_wind_type
        return result


    def _save_to_file(self, datadata, settings):
        print('saving to file: {}'.format(settings['fname']))
        with open(settings['fname'], 'w') as dfile:
            dfile.write(settings['headstring'])
            dfile.write(settings['columnstring'])
            for data in datadata.values():
                # print(data)
                d_date = data['Time']
                d_siteid = data['Site identifier']
                (d_wspd, d_wdir) = polar(data['Value'])
                d_wdir = d_wdir * 57.29577951308232
                if d_wdir < 0: d_wdir = d_wdir + 360.
                dfile.write('{}, {}, {}, {}\n'.format(d_date, d_siteid, d_wspd, d_wdir))
        if self._verbose > 0:
            print('saved')


class MetExtractorPollen(MetExtractor):
    MEASUREMENT_NAME = 'pollen'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollen, self).__init__(out_dir, verbose)
        self._head_string = '{} pollen daily count{} for date range: {} to {}\n'\
            .format(self.MEASUREMENT_NAME, '{}', '{}', '{}')
        self._cols_specific = [self.MEASUREMENT_NAME]
        self._filename = '{}/pollen_{}{}{}.csv'.format(self._out_dir, self.MEASUREMENT_NAME, '{}', '{}')

    @staticmethod
    def get_pollen_species():
        return [s.MEASUREMENT_NAME for s in MetExtractor.all_subclasses(MetExtractorPollen)]


class MetExtractorPollenAlnus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.alnus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenAlnus, self).__init__(out_dir, verbose)


class MetExtractorPollenAmbrosia(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.ambrosia'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenAmbrosia, self).__init__(out_dir, verbose)


class MetExtractorPollenArtemesia(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.artemisia'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenArtemesia, self).__init__(out_dir, verbose)


class MetExtractorPollenBetula(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.betula'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenBetula, self).__init__(out_dir, verbose)


class MetExtractorPollenCorylus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.corylus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenCorylus, self).__init__(out_dir, verbose)


class MetExtractorPollenFraxinus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.fraxinus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenFraxinus, self).__init__(out_dir, verbose)


class MetExtractorPollenPlatanus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.platanus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenPlatanus, self).__init__(out_dir, verbose)


class MetExtractorPollenPoacea(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.poaceae'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenPoacea, self).__init__(out_dir, verbose)


class MetExtractorPollenQuercus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.quercus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenQuercus, self).__init__(out_dir, verbose)


class MetExtractorPollenSalix(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.salix'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenSalix, self).__init__(out_dir, verbose)


class MetExtractorPollenUlmus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.ulmus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenUlmus, self).__init__(out_dir, verbose)


class MetExtractorPollenUrtica(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.urtica'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollenUrtica, self).__init__(out_dir, verbose)


class MetExtractorPollenGroup(object):
    DEFAULT_OUTDIR = 'Pollen_outputs'
    DEFAULT_OUTFILE_SUFFIX = '_pollen'
    DEFAULT_VERBOSE = 0

    def __init__(self, out_dir=DEFAULT_OUTDIR, verbose=DEFAULT_VERBOSE):
        self._extractors = [
            MetExtractorPollenAlnus(out_dir, verbose),
            MetExtractorPollenAmbrosia(out_dir, verbose),
            MetExtractorPollenArtemesia(out_dir, verbose),
            MetExtractorPollenBetula(out_dir, verbose),
            MetExtractorPollenCorylus(out_dir, verbose),
            MetExtractorPollenFraxinus(out_dir, verbose),
            MetExtractorPollenPlatanus(out_dir, verbose),
            MetExtractorPollenPoacea(out_dir, verbose),
            MetExtractorPollenQuercus(out_dir, verbose),
            MetExtractorPollenSalix(out_dir, verbose),
            MetExtractorPollenUlmus(out_dir, verbose),
            MetExtractorPollenUrtica(out_dir, verbose)
        ]
        self._out_dir = out_dir
        self._verbose = verbose

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets=[]):
        result = []
        for pollen_extractor in self._extractors:
            result.append(pollen_extractor.extract_data(
                date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets))
        return result



def create_directory(dir_name):
    try:
        os.makedirs(dir_name)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


if __name__ == '__main__':

    ### general help text
    parser = argparse.ArgumentParser(
        description="*** A script for automated downloading of MET data for a given region and date range. ***")

    AVAILABLE_MEASUREMENTS = MetExtractor.get_top_level_measurements() + MetExtractorPollen.get_pollen_species()


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
    parser.add_argument("--date_range", "-d", dest="date_range", type=str, nargs='+',  help="start and end dates \
                        (array - first two values only). Default: {}".format(str(MetExtractor.DEFAULT_DATE_RANGE)[1:-1].replace(',','')))


    # Latitude / longitude
    parser.add_argument("--latitude_range", "-t", dest="latitude_range", type=int, nargs='+',
                        help="start and end latitude range (array - first two values only). \
                            Default: {}".format(str(MetExtractor.UK_LATITUDES)[1:-1].replace(',','')))
    parser.add_argument("--longitude_range", "-n", dest="longitude_range", type=int, nargs='+',
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
        measurements = args.measurements
        print('Using measurements: {}'.format(measurements))
    else:
        print('No measurements provided, so using default: ', '[' + ", ".join([measure for measure in AVAILABLE_MEASUREMENTS]) + ']')
        measurements = AVAILABLE_MEASUREMENTS

    if args.extra_measurements is True:
        extra_measurements = True
        print('Using extra measurements')
    else:
        print('Not using extra measurements')
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
        if len(args.date_range) >= 2:
            date_range = args.date_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 dates from input --date_range: {}'.format(str(args.date_range)))
        print('Using date range: {}'.format(str(date_range)[1:-1].replace(',','')))
    else:
        print('No date_range provided, so using default: {}'.format(str(MetExtractor.DEFAULT_DATE_RANGE)[1:-1].replace(',','')))
        date_range = MetExtractor.DEFAULT_DATE_RANGE

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


    # Check directory name is valid and create directory
    try:
        print('Creating directory: {}, unless it already exists.'.format(outdir_name))
        create_directory(outdir_name)
    except:
        raise ValueError('Unable to create directory: {}'.format(outdir_name))


    ### Prepare inputs and perform data extraction

    # Create necessary extractors

    # Todo - there should be a factory class for MetExtractor so the classes just get created based on the
    #   measurement input strings.

    if 'rain' in measurements:
        met_extractor_rain = MetExtractorRain(outdir_name, verbose)
        met_extractor_rain.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)
    if 'temp' in measurements:
        met_extractor_temp = MetExtractorTemp(outdir_name, verbose)
        met_extractor_temp.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)
    if 'rel_hum' in measurements:
        met_extractor_relhum = MetExtractorRelativeHumidity(outdir_name, verbose)
        met_extractor_relhum.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)
    if 'pressure' in measurements:
        met_extractor_press = MetExtractorStationPressure(outdir_name, verbose)
        met_extractor_press.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)
    if 'dewpoint' in measurements:
        met_extractor_dewpoint = MetExtractorDewpoint(outdir_name, verbose)
        met_extractor_dewpoint.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)
    if 'wind' in measurements:
        met_extractor_wind = MetExtractorWind(outdir_name, verbose)
        met_extractor_wind.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)
    if 'pollen' in measurements:
        met_extractor_pollen = MetExtractorPollenGroup(outdir_name, verbose)
        met_extractor_pollen.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                          extract_extra_datasets=extra_measurements)
    if 'urtica' in measurements:
        met_extractor_urtica = MetExtractorPollenUrtica(outdir_name, verbose)
        met_extractor_urtica.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
