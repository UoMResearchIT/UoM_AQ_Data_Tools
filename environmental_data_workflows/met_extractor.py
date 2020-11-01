from abc import ABCMeta, abstractmethod
try:
    from medmi_database import Dataset
except:
    print('Warning: Unable to load medmi_database. Not a problem if not downloading MEDMI data.')
from cmath import polar
import json

from .environment_workflow import EnvironmentWorkflow


class MetExtractor(EnvironmentWorkflow):
    __metaclass__ = ABCMeta
    ADDITIONAL_VALS_KEY = 'Additional field values'
    DEFAULT_ADD_EXTRA_MEASUREMENTS = False

    def __init__(self, out_dir=EnvironmentWorkflow.DEFAULT_OUT_DIR, verbose=EnvironmentWorkflow.DEFAULT_VERBOSE):
        super(MetExtractor, self).__init__(out_dir, verbose)
        self._headstring = None
        self._cols_specific = []
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


    def _get_extraction_dict(self):
        result = {
            'Time range': self.date_range,
            'Source reference': self.SOURCE_REFERENCE,
            'Latitude range': self.latitude_range,
            'Longitude range': self.longitude_range
        }

        return result

    def _get_settings(self, outfile_suffix=EnvironmentWorkflow.DEFAULT_OUT_FILE_SUFFIX):
        str_extra_datasets = '-'.join(self._extra_datasets).replace(' ', '-')
        if len(self._extra_datasets) > 0:
            filename = self._file_out.format('_extras-{}'.format(str_extra_datasets), outfile_suffix)
            headstring = self._head_string.format(' and extra datasets: {}'.format(str(str_extra_datasets)),
                                                   self.date_range[0], self.date_range[1])
        else:
            filename = self._file_out.format('', outfile_suffix)
            headstring = self._head_string.format('', self.date_range[0], self.date_range[1])

        result = {
            'fname': filename,
            'headstring': headstring,
            'columnstring': ','.join(self.get_all_cols()) + '\n'
        }
        return result

    def get_all_cols(self):
        # Met extractor has extra_datasets columns, in addition to the usual (base and specific) columns.
        return super(MetExtractor, self).get_all_cols() + self._extra_datasets


    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets=False):
        self.set_region(latitude_range, longitude_range)
        self.date_range = date_range

        if extract_extra_datasets:
            self._extra_datasets = self.ALLOWED_EXTRA_DATASETS
        else:
            self._extra_datasets = []

        extraction_dict = self._get_extraction_dict()
        settings = self._get_settings(outfile_suffix)
        print('extracting {}'.format(settings['headstring'], date_range[0], date_range[1]))
        return self._extract_data(extraction_dict, settings)


    def _extract_data(self, extraction_dict, settings, save_to_file=True):
        if self.verbose >= 1:
            print('extracting data for {}'.format(self.MEASUREMENT_NAME))
        if self.verbose > 1:
            print('using extraction dict: {}'.format(json.dumps(extraction_dict)))
            print('using settings: {}'.format(json.dumps(settings)))
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
        if self.verbose > 0:
            print('saved')


class MetExtractorRain(MetExtractor):
    SOURCE_REFERENCE = 'midas.rain_drnl_ob.prcp_amt 1'
    MEASUREMENT_NAME = 'rain'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorRain, self).__init__(out_dir, verbose)
        self._head_string = 'Rain gauge daily data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorRain.MEASUREMENT_NAME]
        self._file_out = '{}/rain{}{}.csv'.format(self.out_dir, '{}', '{}')


class MetExtractorTemperature(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.air_temperature'
    MEASUREMENT_NAME = 'temperature'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'pressure', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorTemperature, self).__init__(out_dir, verbose)
        self._head_string = 'Temperature data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorTemperature.MEASUREMENT_NAME]
        self._file_out = '{}/temp{}{}.csv'.format(self.out_dir, '{}', '{}')


class MetExtractorRelativeHumidity(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.rltv_hum'
    MEASUREMENT_NAME = 'rel_hum'
    ALLOWED_EXTRA_DATASETS = ['temperature', 'pressure', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorRelativeHumidity, self).__init__(out_dir, verbose)
        self._head_string = 'Relative Humidity data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorRelativeHumidity.MEASUREMENT_NAME]
        self._file_out = '{}/rel_hum{}{}.csv'.format(self.out_dir, '{}', '{}')


class MetExtractorStationPressure(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.stn_pres'
    MEASUREMENT_NAME = 'pressure'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'temperature', 'dewpoint']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorStationPressure, self).__init__(out_dir, verbose)
        self._head_string = 'Station Pressure data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorStationPressure.MEASUREMENT_NAME]
        self._file_out = '{}/pressure{}{}.csv'.format(self.out_dir, '{}', '{}')


class MetExtractorDewpoint(MetExtractor):
    SOURCE_REFERENCE = 'midas.weather_hrly_ob.dewpoint'
    MEASUREMENT_NAME = 'dewpoint'
    ALLOWED_EXTRA_DATASETS = ['rel_hum', 'pressure', 'temperature']

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorDewpoint, self).__init__(out_dir, verbose)
        self._head_string = 'Dewpoint hourly data{} for date range: {} to {}\n'
        self._cols_specific = [MetExtractorDewpoint.MEASUREMENT_NAME]
        self._file_out = '{}/dewpoint{}{}.csv'.format(self.out_dir, '{}', '{}')


class MetExtractorWind(MetExtractor):
    SOURCE_REFERENCE = 'midas.wind_mean_ob.mean_wind_speed 1'
    MEASUREMENT_NAME = 'wind'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorWind, self).__init__(out_dir, verbose)
        self._complex_wind_type = True
        self._head_string = 'Wind speed/direction hourly data{} for date range: {} to {}\n'
        self._cols_specific = ['windspeed','winddir']
        self._file_out = '{}/wind{}{}.csv'.format(out_dir, '{}', '{}')

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
        if self.verbose > 0:
            print('saved')


class MetExtractorPollen(MetExtractor):
    MEASUREMENT_NAME = 'pollen'
    ALLOWED_EXTRA_DATASETS = []

    def __init__(self, out_dir=MetExtractor.DEFAULT_OUT_DIR, verbose=MetExtractor.DEFAULT_VERBOSE):
        super(MetExtractorPollen, self).__init__(out_dir, verbose)
        self._head_string = '{} pollen daily count{} for date range: {} to {}\n'\
            .format(self.MEASUREMENT_NAME, '{}', '{}', '{}')
        self._cols_specific = [self.MEASUREMENT_NAME]
        self._file_out = '{}/pollen_{}{}{}.csv'.format(self.out_dir, self.MEASUREMENT_NAME, '{}', '{}')
        if self.MEASUREMENT_NAME == 'pollen':
            self._extractor = MetExtractorPollenGroup(out_dir, verbose)
        else:
            self._extractor = super(MetExtractorPollen, self)

    @staticmethod
    def get_pollen_species():
        return [s.MEASUREMENT_NAME for s in MetExtractor.all_subclasses(MetExtractorPollen)]

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets=[]):
        return self._extractor.extract_data(date_range, latitude_range, longitude_range, outfile_suffix,
                                            extract_extra_datasets)


class MetExtractorPollenAlnus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.alnus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenAmbrosia(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.ambrosia'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenArtemesia(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.artemisia'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]

class MetExtractorPollenBetula(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.betula'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenCorylus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.corylus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenFraxinus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.fraxinus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenPlatanus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.platanus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenPoacea(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.poaceae'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenQuercus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.quercus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenSalix(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.salix'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenUlmus(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.ulmus'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenUrtica(MetExtractorPollen):
    SOURCE_REFERENCE = 'midas.pollen_drnl_ob.urtica'
    MEASUREMENT_NAME = SOURCE_REFERENCE.split('.')[-1]


class MetExtractorPollenGroup(object):
    DEFAULT_OUTDIR = 'Pollen_outputs'
    DEFAULT_OUTFILE_SUFFIX = '_pollen'
    DEFAULT_VERBOSE = 0

    def __init__(self, out_dir=DEFAULT_OUTDIR, verbose=DEFAULT_VERBOSE):
        self._extractors = []
        for species in MetExtractorPollen.get_pollen_species():
            self._extractors.append(MetExtractor.get_class_from_measurement_name(species)(out_dir, verbose))
        self.out_dir = out_dir
        self.verbose = verbose

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets=[]):
        result = []
        for pollen_extractor in self._extractors:
            result.append(pollen_extractor.extract_data(
                date_range, latitude_range, longitude_range, outfile_suffix, extract_extra_datasets))
        return result