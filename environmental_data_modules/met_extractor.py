from abc import ABCMeta, abstractmethod
try:
    from medmi_database import Dataset
except:
    pass
from cmath import polar
from datetime import datetime
import json

from environmental_data_modules import Extractor, MetModule, DateRangeProcessor, RegionRectProcessor


class MetExtractor(Extractor, MetModule, DateRangeProcessor, RegionRectProcessor):
    """
        Abstract class, parent of classes used for extracting data from the MEDMI server.
    """
    __metaclass__ = ABCMeta

    # Define 'absolute' constants
    ADDITIONAL_VALS_KEY = 'Additional field values'
    BASE_FILE_OUT = '{}/Met_extracted_{}{}{}.csv'
    MEDMI_DATE_FORMAT = '%Y-%m-%d %H'

    # Define default constants
    DEFAULT_ADD_EXTRA_MEASUREMENTS = False

    def __init__(self, out_dir=Extractor.DEFAULT_OUT_DIR, verbose=Extractor.DEFAULT_VERBOSE):
        """ Initialise instance of the MetExtractor class.
            Initialises the private class variables

            Args:
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of MetExtractor

        """
        super(MetExtractor, self).__init__(out_dir, verbose)
        MetModule.__init__(self)
        DateRangeProcessor.__init__(self)
        RegionRectProcessor.__init__(self)
        self._headstring = None
        self._base_file_out = MetExtractor.BASE_FILE_OUT
        self._temp_file_out = ''
        self._extra_datasets = []

    @staticmethod
    def get_source_ref_from_name(name):
        for sub_class in MetExtractor.all_subclasses(MetExtractor):
            if hasattr(sub_class, 'MEASUREMENT_NAME') and sub_class.MEASUREMENT_NAME == name:
                return sub_class.SOURCE_REFERENCE
        raise ValueError('Measurement does not exist: {}'.format(name))

    @staticmethod
    def get_class_from_measurement_name(name):
        for sub_class in MetExtractor.all_subclasses(MetExtractor):
            if hasattr(sub_class, 'MEASUREMENT_NAME') and sub_class.MEASUREMENT_NAME == name:
                return sub_class
        raise ValueError('Measurement does not exist: {}'. format(name))

    @staticmethod
    def get_available_measurements():
        return [s.MEASUREMENT_NAME for s in MetExtractor.all_subclasses(MetExtractor)]

    def get_all_column_headers(self):
        return super(MetExtractor, self).get_all_column_headers() + self._extra_datasets

    def _get_extraction_dict(self):
        result = {
            'Time range': [self.date_range[0].strftime(MetExtractor.MEDMI_DATE_FORMAT),
                           self.date_range[1].strftime(MetExtractor.MEDMI_DATE_FORMAT)],
            'Source reference': self.SOURCE_REFERENCE,
            'Latitude range': self.latitude_range,
            'Longitude range': self.longitude_range
        }

        return result

    def extract_data(self, date_range=None,
                     latitude_range=RegionRectProcessor.UK_LATITUDES,
                     longitude_range=RegionRectProcessor.UK_LONGITUDES,
                     outfile_suffix=Extractor.DEFAULT_OUT_FILE_SUFFIX,
                     extract_extra_datasets=DEFAULT_ADD_EXTRA_MEASUREMENTS):

        """ Extract the date from MEDMI server, based on parameters

            Args:
                date_range:         (list of 2 datetime) The date range of interest
                latitude_range:     (list of 2 floats) The latitude range of interest
                longitude_range:    (list of 2 floats) The longitude range of interest
                outfile_suffix:     (string) The suffix to appended to the end of output file names.
                extract_extra_datasets: (boolean). Whether to extract the available extra data measurements.

            Returns:
                Extracted data (pandas dataframe)

        """

        self.set_region(latitude_range, longitude_range)
        self._outfile_suffix = outfile_suffix

        if date_range is not None:
            self.date_range = [datetime.strptime(date_range[0], DateRangeProcessor.INPUT_DATE_FORMAT),
                               datetime.strptime(date_range[1], DateRangeProcessor.INPUT_DATE_FORMAT)]
        else:
            self.date_range = [self.get_available_start(), self.get_available_end()]

        if extract_extra_datasets:
            self._extra_datasets = self.ALLOWED_EXTRA_DATASETS
        else:
            self._extra_datasets = []

        extraction_dict = self._get_extraction_dict()

        if len(self._extra_datasets) > 0:
            str_extra_datasets = '-'.join(self._extra_datasets).replace(' ', '-')
            self.file_out = self._temp_file_out.format('_extras-{}'.format(str_extra_datasets),
                                                       self.outfile_suffix_string)
            self._headstring = self._head_string.format(' and extra datasets: {}'.format(str(str_extra_datasets)),
                                                        self.start, self.end)
        else:
            self.file_out = self._temp_file_out.format('', self.outfile_suffix_string)
            self._headstring = self._head_string.format('', self.start, self.end)

        print('extracting {}'.format(self._headstring, self.start, self.end))
        return self._extract_data(extraction_dict)

    def _extract_data(self, extraction_dict, save_to_csv=True):
        if self.verbose >= 1:
            print('extracting data for {}'.format(self.MEASUREMENT_NAME))
        if self.verbose > 1:
            print('using extraction dict: {}'.format(json.dumps(extraction_dict, default=str)))
        datadata = self._perform_extraction(extraction_dict)
        if save_to_csv:
            self._save_to_csv(datadata)
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

    def _save_to_csv(self, data_result):
        print('saving to file: {}'.format(self.file_out))

        with open(self.file_out, 'w') as dfile:
            dfile.write(self._headstring)
            dfile.write(','.join(self.get_all_column_headers()) + '\n')

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
        if self.verbose > 0:
            print('saved')


class MetExtractorRain(MetExtractor):
    """
        Class used for extracting rain data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.rain_drnl_ob.prcp_amt 1'
    MEASUREMENT_NAME = 'rain'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorRain, self).__init__(out_dir, verbose)
        self._head_string = 'Rain gauge daily data{} for date range: {} to {}\n'
        self._columns_specific = [MetExtractorRain.MEASUREMENT_NAME]
        self._temp_file_out = self.base_file_out.format(self.out_dir, 'rain', '{}', '{}')


class MetExtractorTemperature(MetExtractor):
    """
        Class used for extracting temperature data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.air_temperature'
    MEASUREMENT_NAME = 'temperature'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'pressure', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorTemperature, self).__init__(out_dir, verbose)
        self._head_string = 'Temperature data{} for date range: {} to {}\n'
        self._columns_specific = [MetExtractorTemperature.MEASUREMENT_NAME]
        self._temp_file_out = self.base_file_out.format(self.out_dir, 'temp', '{}', '{}')


class MetExtractorRelativeHumidity(MetExtractor):
    """
        Class used for extracting relative humidity data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.rltv_hum'
    MEASUREMENT_NAME = 'rel_hum'
    ALLOWED_EXTRA_DATASETS = ['temperature', 'pressure', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorRelativeHumidity, self).__init__(out_dir, verbose)
        self._head_string = 'Relative Humidity data{} for date range: {} to {}\n'
        self._columns_specific = [MetExtractorRelativeHumidity.MEASUREMENT_NAME]
        self._temp_file_out = self.base_file_out.format(self.out_dir, 'rel_hum', '{}', '{}')


class MetExtractorStationPressure(MetExtractor):
    """
        Class used for extracting station pressure data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.stn_pres'
    MEASUREMENT_NAME = 'pressure'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'temperature', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorStationPressure, self).__init__(out_dir, verbose)
        self._head_string = 'Station Pressure data{} for date range: {} to {}\n'
        self._columns_specific = [MetExtractorStationPressure.MEASUREMENT_NAME]
        self._temp_file_out = self.base_file_out.format(self.out_dir, 'pressure', '{}', '{}')


class MetExtractorDewpoint(MetExtractor):
    """
        Class used for extracting dewpoint data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.dewpoint'
    MEASUREMENT_NAME = 'dewpoint'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'pressure', 'temperature']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorDewpoint, self).__init__(out_dir, verbose)
        self._head_string = 'Dewpoint hourly data{} for date range: {} to {}\n'
        self._columns_specific = [MetExtractorDewpoint.MEASUREMENT_NAME]
        self._temp_file_out = self.base_file_out.format(self.out_dir, 'dewpoint', '{}', '{}')


class MetExtractorWind(MetExtractor):
    """
        Class used for extracting wind data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.wind_mean_ob.mean_wind_speed 1'
    MEASUREMENT_NAME = 'wind'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorWind, self).__init__(out_dir, verbose)
        self._complex_wind_type = True
        self._head_string = 'Wind speed/direction hourly data{} for date range: {} to {}\n'
        self._columns_specific = ['windspeed','winddir']
        self._temp_file_out = self.base_file_out.format(self.out_dir, 'wind', '{}', '{}')

    def _get_extraction_dict(self):
        result = super(MetExtractorWind, self)._get_extraction_dict()
        result['Complex wind type'] = self._complex_wind_type
        return result

    def _save_to_csv(self, datadata):
        print('saving to file: {}'.format(self.file_out))
        with open(self.file_out, 'w') as dfile:
            dfile.write(self._headstring)
            dfile.write(','.join(self.get_all_column_headers()) + '\n')
            for data in datadata.values():
                # print(data)
                d_date = data['Time']
                d_siteid = data['Site identifier']
                (d_wspd, d_wdir) = polar(data['Value'])
                d_wdir = d_wdir * 57.29577951308232
                if d_wdir < 0: d_wdir = d_wdir + 360.
                dfile.write('{}, {}, {}, {}\n'.format(d_date, d_siteid, d_wspd, d_wdir))
        if self.verbose > 0:
            print('saved')


class MetExtractorPollen(MetExtractor):
    """
        Class used for extracting all pollen species data from the MEDMI server.
        Parent class of individual pollen species extractors (e.g. urtica)
    """
    MEASUREMENT_NAME = 'pollen'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollen, self).__init__(out_dir, verbose)
        self._head_string = '{} pollen daily count{} for date range: {} to {}\n'\
            .format(self.MEASUREMENT_NAME, '{}', '{}', '{}')
        self._columns_specific = [self.MEASUREMENT_NAME]
        self._temp_file_out = self.base_file_out.format(self.out_dir, 'pollen-{}'.format(self.MEASUREMENT_NAME), '{}', '{}')
        if self.MEASUREMENT_NAME == 'pollen':
            self._extractor = MetExtractorPollenGroup(out_dir, verbose)
        else:
            self._extractor = super(MetExtractorPollen, self)

    @staticmethod
    def get_pollen_species():
        return [s.MEASUREMENT_NAME for s in MetExtractor.all_subclasses(MetExtractorPollen)]

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets):
        return self._extractor.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                            extract_extra_datasets)


class MetExtractorPollenAlnus(MetExtractorPollen):
    """
        Class used for extracting alnus pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.alnus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenAmbrosia(MetExtractorPollen):
    """
        Class used for extracting ambrosia pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.ambrosia'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenArtemesia(MetExtractorPollen):
    """
        Class used for extracting artemesia pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.artemisia'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenBetula(MetExtractorPollen):
    """
        Class used for extracting betula pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.betula'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenCorylus(MetExtractorPollen):
    """
        Class used for extracting corylus pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.corylus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenFraxinus(MetExtractorPollen):
    """
        Class used for extracting fraxinus pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.fraxinus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenPlatanus(MetExtractorPollen):
    """
        Class used for extracting platanus pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.platanus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenPoacea(MetExtractorPollen):
    """
        Class used for extracting poacea pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.poaceae'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenQuercus(MetExtractorPollen):
    """
        Class used for extracting quercus pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.quercus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenSalix(MetExtractorPollen):
    """
        Class used for extracting salix pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.salix'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenUlmus(MetExtractorPollen):
    """
        Class used for extracting ulmus pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.ulmus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenUrtica(MetExtractorPollen):
    """
        Class used for extracting urtica pollen data from the MEDMI server.
    """
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.urtica'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenGroup(object):
    """
        Class used to contain a collection of pollen data extractors for extracting from the MEDMI server.
    """
    DEFAULT_OUTDIR = 'Pollen_outputs'
    DEFAULT_OUTFILE_SUFFIX = '_pollen'
    DEFAULT_VERBOSE = 0

    def __init__(self, out_dir=DEFAULT_OUTDIR, verbose=DEFAULT_VERBOSE):
        self._extractors = []
        for species in MetExtractorPollen.get_pollen_species():
            self._extractors.append(MetExtractor.get_class_from_measurement_name(species)(out_dir, verbose))
        self.out_dir = out_dir
        self.verbose = verbose

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets):
        result = []
        for pollen_extractor in self._extractors:
            result.append(pollen_extractor.extract_data(
                date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets))
        return result
