from abc import ABCMeta, abstractmethod
from medmi_database import Dataset
from cmath import polar
import argparse
import json

import os, errno


class MetExtractor(object):
    __metaclass__ = ABCMeta
    VERBOSE = 0
    DEFAULT_OUT_FILE_SUFFIX = ''
    DEFAULT_OUT_DIR = 'met_extracted_data'
    DEFAULT_DATE_RANGE = ['2016-1-1 0', '2019-12-31 23']
    UK_LATITUDES = [48, 60]
    UK_LONGITUDES = [-11, 3]
    DEFAULT_VERBOSE = 0
    DEFAULT_PROCESS_SPATIAL = 'sp_idw_mean'
    DEFAULT_PROCESS_TEMPORAL = 'tp_mean'
    DEFAULT_PROCESS_RADIUS = 50000
    DEFAULT_PROCESS_PERIOD = 1
    DEFAULT_COLS_BASE = 'date,siteID,{}\n'
    ADDITIONAL_VALS_KEY = 'Additional field values'



    def __init__(self, dir_name=DEFAULT_OUT_DIR, verbose=DEFAULT_VERBOSE):
        self._out_dir = dir_name
        self._spatial_process = None
        self._temporal_process = None
        self._cols_base = MetExtractor.DEFAULT_COLS_BASE
        self._latitude_range = MetExtractor.UK_LATITUDES
        self._longitude_range = MetExtractor.UK_LONGITUDES
        self._longitude = None
        self._latitude = None
        self._date_range = None
        self._source_reference = None
        self._filename = None
        self._headstring = None
        self._cols_specific = None
        self._verbose=verbose

    @abstractmethod
    def get_source_reference(self):
        raise NotImplementedError("Must override get_source_reference")

    @abstractmethod
    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        raise NotImplementedError("Must override extract_data")

    @abstractmethod
    def extract_data_from_dict(self, extraction_dict, settings, save_to_file):
        raise NotImplementedError("Must override extract_data_from_dict")

    def _apply_spatial_process(self, dataset):
        dataset.process(self._spatial_process)

    def add_outfile_suffix(self, outfile_suffix):
        return self._filename.format(outfile_suffix)

    def get_settings(self):
        result = {
            'fname': self._filename,
            'headstring': self._head_string,
            'columnstring': self._cols_base.format(self._cols_specific)
        }
        return result

    def _apply_temporal_process(self, dataset):
        dataset.process(self._temporal_process)

    def set_spatial_process(self, method=DEFAULT_PROCESS_SPATIAL, radius=DEFAULT_PROCESS_RADIUS):
        self._spatial_process = {'Method': method, 'Radius': radius}

    def set_temporal_process(self, method=DEFAULT_PROCESS_TEMPORAL, period=DEFAULT_PROCESS_PERIOD):
        self._temporal_process = {'Method': method, 'Period': period}

    def set_date_range(self, date_range):
        self._date_range = date_range

    def get_date_range(self):
        return self._date_range

    def get_out_dir(self):
        return self._out_dir

    def set_region(self, latitude_range, longitude_range):
        self._latitude_range = latitude_range
        self._longitude_range = longitude_range

        # Can only do regions or sites
        self._latitude = None
        self._longitude = None

    def set_site(self, latitude, longitude):
        self._latitude = latitude
        self._longitude = longitude

        # Can only do regions or sites
        self._longitude_range = None
        self._latitude_range = None

    def _get_extraction_dict(self):
        result = {
            'Time range': self._date_range,
            'Source reference': self.get_source_reference()
        }
        if self._latitude is not None:
            result['Latitude'] = self._latitude
            result['Longitude'] = self._longitude
        elif self._latitude_range is not None:
            result['Latitude range'] =  self._latitude_range
            result['Longitude range'] = self._longitude_range
        return result

    def extract_data_from_dict(self, extraction_dict, settings, save_to_file=True):
        return self._extract_data(extraction_dict, settings, save_to_file=save_to_file)


    def _extract_data(self, dict, settings, save_to_file=True, extra_datasets=[]):
        if self.VERBOSE == 1:
            print('extracting data')
        elif self.VERBOSE > 1:
            print('extracting data using dict: {}\n and settings: {}'.format(json.dumps(dict, default=str), format(json.dumps(settings))))

        datadata = self._perform_extraction(dict)

        for ds in extra_datasets:
            datadata.add(ds)

        if save_to_file:
            self.save_to_file(datadata, settings)

        return datadata

    def save_to_file(self, datadata, settings):
        print('saving to file: {}'.format(settings['fname']))
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
        print('saved')

    def _perform_extraction(self, dict):
        datadata = Dataset(dict)
        if self._spatial_process is not None:
            if self.VERBOSE > 0:
                print('applying spatial processing: {}'.format(self._spatial_process))
            self._apply_spatial_process(datadata)
        if self._temporal_process is not None:
            if self.VERBOSE > 0:
                print('applying temporal processing: {}'.format(self._temporal_process))
            self._apply_temporal_process(datadata)
        if self._temporal_process is None and self._spatial_process is None:
            if self.VERBOSE > 0:
                print('applying default')
            datadata.default()
        return datadata



class MetExtractorRain(MetExtractor):
    def __init__(self, outfile_dir=MetExtractor.DEFAULT_OUT_DIR):
        super(MetExtractorRain, self).__init__(outfile_dir)
        self._source_reference = 'midas.rain_drnl_ob.prcp_amt 1'
        self._head_string = 'Rain gauge daily data for date range: {} to {}\n'
        self._cols_specific = 'rain'
        self._filename = '{}/rain{}.csv'.format(self._out_dir, '{}')

    def get_source_reference(self):
        return self._source_reference

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        super(MetExtractorRain, self).set_region(latitude_range, longitude_range)
        super(MetExtractorRain, self).set_date_range(date_range)

        print('extracting {}'.format(self._head_string.format(date_range[0], date_range[1])))
        extraction_dict = super(MetExtractorRain, self)._get_extraction_dict()
        settings = {'fname': self._filename.format(outfile_suffix),
                    'headstring': self._head_string.format(date_range[0], date_range[1]),
                    'columnstring': self._cols_base.format(self._cols_specific)}
        return self._extract_data(extraction_dict, settings)


class MetExtractorTemp(MetExtractor):
    DEFAULT_PROCESS_SPATIAL_TEMPERATURE = 'sp_altadj_idw_mean'

    def __init__(self, outfile_dir=MetExtractor.DEFAULT_OUT_DIR):
        super(MetExtractorTemp, self).__init__(outfile_dir)
        self._source_reference = 'midas.weather_hrly_ob.air_temperature'
        self._head_string = 'Temperature, Relative Humidity, Station Pressure, and Wet Bulb Temperature hourly \
                                    data for date range: {} to {}\n'
        self._cols_specific = 'temperature'
        self._spatial_temperature_process = None
        self._filename = '{}/temp_rh_press_wbulb{}.csv'.format(self._out_dir, '{}')

    def _apply_spatial_process(self, dataset):
        dataset.process(self._spatial_temperature_process)

    def get_source_reference(self):
        return self._source_reference

    def set_spatial_temperature_process(self, method=DEFAULT_PROCESS_SPATIAL_TEMPERATURE,
                                        radius=MetExtractor.DEFAULT_PROCESS_RADIUS):
        self._spatial_temperature_process = {'Method': method, 'Period': radius}

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix=''):
        super(MetExtractorTemp, self).set_region(latitude_range, longitude_range)
        super(MetExtractorTemp, self).set_date_range(date_range)
        print('extracting {}'.format(self._head_string.format(date_range[0], date_range[1])))
        extraction_dict = super(MetExtractorTemp, self)._get_extraction_dict()

        settings = {'fname': self._filename.format(outfile_suffix),
                    'headstring': self._head_string.format(date_range[0], date_range[1]),
                    'columnstring': self._cols_base.format(self._cols_specific)}
        return self._extract_data(extraction_dict, settings)


class MetExtractorRelativeHumidity(MetExtractor):
    def __init__(self, outfile_dir=MetExtractor.DEFAULT_OUT_DIR):
        super(MetExtractorRelativeHumidity, self).__init__(outfile_dir)
        self._source_reference = 'midas.weather_hrly_ob.rltv_hum'
        self._head_string = 'Relative Humidity  data for date range: {} to {}\n'
        self._cols_specific = 'rh'
        self._filename = '{}/rh{}.csv'.format(self._out_dir, '{}')

    def get_source_reference(self):
        return self._source_reference

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix=''):
        super(MetExtractorRelativeHumidity, self).set_region(latitude_range, longitude_range)
        super(MetExtractorRelativeHumidity, self).set_date_range(date_range)
        print('extracting {}'.format(self._head_string.format(date_range[0], date_range[1])))
        extraction_dict = super(MetExtractorRelativeHumidity, self)._get_extraction_dict()

        settings = {'fname': self._filename.format(outfile_suffix),
                    'headstring': self._head_string.format(date_range[0], date_range[1]),
                    'columnstring': self._cols_base.format(self._cols_specific)}
        return self._extract_data(extraction_dict, settings)


class MetExtractorStationPressure(MetExtractor):
    def __init__(self, outfile_dir=MetExtractor.DEFAULT_OUT_DIR):
        super(MetExtractorStationPressure, self).__init__(outfile_dir)
        self._source_reference = 'midas.weather_hrly_ob.stn_pres'
        self._head_string = 'Station Pressure data for date range: {} to {}\n'
        self._cols_specific = 'pressure'
        self._filename = '{}/pressure{}.csv'.format(self._out_dir, '{}')

    def get_source_reference(self):
        return self._source_reference

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix=''):
        super(MetExtractorStationPressure, self).set_region(latitude_range, longitude_range)
        super(MetExtractorStationPressure, self).set_date_range(date_range)
        print('extracting {}'.format(self._head_string.format(date_range[0], date_range[1])))
        extraction_dict = super(MetExtractorStationPressure, self)._get_extraction_dict()

        settings = {'fname': self._filename.format(outfile_suffix),
                    'headstring': self._head_string.format(date_range[0], date_range[1]),
                    'columnstring': self._cols_base.format(self._cols_specific)}
        return self._extract_data(extraction_dict, settings)


class MetExtractorDewpoint(MetExtractor):

    def __init__(self, outfile_dir=MetExtractor.DEFAULT_OUT_DIR):
        super(MetExtractorDewpoint, self).__init__(outfile_dir)
        self._source_reference = 'midas.weather_hrly_ob.dewpoint'
        self._head_string = 'Wet Bulb Temperature hourly data for date range: {} to {}\n'
        self._cols_specific = 'wbtemp'
        self._spatial_temperature_process = None
        self._filename = '{}/wet_bulb{}.csv'.format(self._out_dir, '{}')

    def get_source_reference(self):
        return self._source_reference

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix=''):
        super(MetExtractorDewpoint, self).set_region(latitude_range, longitude_range)
        super(MetExtractorDewpoint, self).set_date_range(date_range)
        print('extracting {}'.format(self._head_string.format(date_range[0], date_range[1])))
        extraction_dict = super(MetExtractorDewpoint, self)._get_extraction_dict()

        settings = {'fname': self._filename.format(outfile_suffix),
                    'headstring': self._head_string.format(date_range[0], date_range[1]),
                    'columnstring': self._cols_base.format(self._cols_specific)}
        return self._extract_data(extraction_dict, settings)


class MetExtractorWind(MetExtractor):

    def __init__(self, outfile_dir=MetExtractor.DEFAULT_OUT_DIR):
        super(MetExtractorWind, self).__init__(outfile_dir)
        self._source_reference = 'midas.wind_mean_ob.mean_wind_speed 1'
        self._complex_wind_type = True
        self._head_string = 'Wind speed and direction hourly data for date range: {} to {}\n'
        self._cols_specific = 'windspeed,winddir'
        self._filename = '{}/wind{}.csv'.format(outfile_dir, '{}')

    def get_source_reference(self):
        return self._source_reference

    def _get_extraction_dict(self):
        result = super(MetExtractorWind, self)._get_extraction_dict()
        result['Complex wind type'] = self._complex_wind_type
        return result

    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        super(MetExtractorWind, self).set_region(latitude_range, longitude_range)
        super(MetExtractorWind, self).set_date_range(date_range)
        print('extracting {}'.format(self._head_string.format(date_range[0], date_range[1])))
        extraction_dict = self._get_extraction_dict()
        settings = {    'fname': self._filename.format(outfile_suffix),
                        'headstring': self._head_string.format(date_range[0], date_range[1]),
                        'columnstring': self._cols_base.format(self._cols_specific)}

        return self._extract_data(extraction_dict, settings)

    def extract_data_from_dict(self, extraction_dict, settings, save_to_file=True):
        extraction_dict['Complex wind type'] = self._complex_wind_type
        return self._extract_data(extraction_dict, settings, save_to_file=save_to_file)

    # Overrides base class method
    def _extract_data(self, dict, settings, save_to_file=True):
        print('extracting using dict: {}\n and settings: {}'.format(json.dumps(dict), format(json.dumps(settings))))
        datadata = self._perform_extraction(dict)
        if save_to_file:
            self.save_to_file(datadata, settings)
        return datadata

    def save_to_file(self, datadata, settings):
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
        print('saved')


class MetExtractorPollen(MetExtractor):
    pollen_sites = {
        'BELFAST': (54.58537, -5.93784, 15.0),
        'BEVERLEY': (53.84243, -0.4252, 25.0),
        'CARDIFF': (51.49536, -3.21084, 25.0),
        'CHESTER': (53.19835, -2.89674, 28.0),
        'ESKDALEMUIR': (55.31184, -3.20545, 236.0),
        'EXETER': (50.73599, -3.53101, 10.0),
        'IPSWICH': (52.05638, 1.1997, 41.0),
        'LEICESTER': (52.62155, -1.12227, 91.0),
        'LONDON': (51.51022, -0.11533, 45.0),
        'MYLNEFIELD': (56.45699, -3.07182, 31.0),
        'PLYMOUTH': (50.37461, -4.13725, 45.0),
        'WIGHT': (50.71052, -1.29944, 32.0),
        'WORCESTER': (52.19716, -2.24165, 40.0),
        'YORK': (53.94819, -1.05194, 25.0)}
    pollen_species = ['alnus', 'ambrosia', 'artemisia', 'betula', 'corylus', 'fraxinus', 'platanus',
                           'poaceae', 'quercus', 'salix', 'ulmus', 'urtica']

    def __init__(self, outfile_dir=MetExtractor.DEFAULT_OUT_DIR):
        super(MetExtractorPollen, self).__init__(outfile_dir)

        self._source_reference_base = 'midas.pollen_drnl_ob.{}'
        self._head_string = '{} pollen daily count for date range: {} to {}\n'
        self._cols_specific = str(MetExtractorPollen.pollen_species)[1:-1]
        self._filename = '{}/pollen_{}{}.csv'.format(self._out_dir, '{}', '{}')

    def get_source_reference(self):
        return self._source_reference_base

    def add_outfile_suffix(self, outfile_suffix):
        return self._filename.format('{}', outfile_suffix)


    def extract_data(self, date_range, latitude_range, longitude_range, outfile_suffix):
        result = []
        super(MetExtractorPollen, self).set_region(latitude_range, longitude_range)
        super(MetExtractorPollen, self).set_date_range(date_range)
        fname = self._filename.format('{}', outfile_suffix)
        for pspc in MetExtractorPollen.pollen_species:
            print('working on pollen species: {}'.format(pspc))
            fname = fname.format(pspc)
            extraction_dict = super(MetExtractorPollen, self)._get_extraction_dict()
            extraction_dict['Source reference'] = self._source_reference_base.format(pspc)
            settings = {'fname': fname,
                        'headstring': self._head_string.format(pspc, self._date_range[0], self._date_range[1]),
                        'columnstring': self._cols_base.format(pspc)}
            result.append(super(MetExtractorPollen, self)._extract_data(extraction_dict, settings))
        return result

    def extract_data_from_dict(self, extraction_dict, settings, save_to_file=True):
        result = []
        for pspc in MetExtractorPollen.pollen_species:
            print('working on pollen species: {}'.format(pspc))
            fname = settings['fname'].format(pspc)
            extraction_dict['Source reference'] = self._source_reference_base.format(pspc)
            settings = {'fname': fname,
                        'headstring': settings['headstring'].format(pspc, self._date_range[0], self._date_range[1]),
                        'columnstring': settings['columnstring'].format(pspc)}
            result.append(super(MetExtractorPollen, self)._extract_data(extraction_dict, settings, save_to_file=save_to_file))
        return result

class MetExtractorGroup(object):
    DEFAULT_METHOD = 'MEDMI TS'
    DEFAULT_PERIOD = 1
    DEFAULT_OUTDIR = 'GROUP_DATA'
    DEFAULT_OUTFILE_SUFFIX = '_site_group'
    DEFAULT_WRITE_INTERIM_OUTFILES = False
    DEFAULT_VERBOSE = 0

    def __init__(self, extractors, out_dir=DEFAULT_OUTDIR, verbose=DEFAULT_VERBOSE):
        self._extractors = extractors
        self._cols_base = 'date{}'
        self._method = MetExtractorGroup.DEFAULT_METHOD
        self._period = MetExtractorGroup.DEFAULT_PERIOD
        self._out_dir = out_dir
        self._filename_base = '{}/site_{}_{}-{}{}.txt'
        self._verbose = VERBOSE

    def get_date_range(self):
        try:
            return self._extractors[0].get_date_range()
        except:
            return None

    def set_all_date_ranges(self, date_range):
        for extractor in self._extractors:
            extractor.set_date_range(date_range)

    def set_all_regions(self, latitude_range, longitude_range):
        for extractor in self._extractors:
            extractor.set_region(latitude_range, longitude_range)

    def set_all_sites(self, latitude, longitude):
        for extractor in self._extractors:
            extractor.set_site(latitude, longitude)

    def set_spatial_processes(self, method=MetExtractor.DEFAULT_PROCESS_SPATIAL,
                            radius=MetExtractor.DEFAULT_PROCESS_RADIUS):
        for extractor in self._extractors:
            extractor.set_spatial_process(method, radius)

    def set_spatial_temperature_process(self, method=MetExtractorTemp.DEFAULT_PROCESS_SPATIAL_TEMPERATURE,
                            radius=MetExtractor.DEFAULT_PROCESS_RADIUS):
        for extractor in self._extractors:
            if type(extractor).__name__ == MetExtractorTemp.__name__:
                extractor.set_spatial_temperature_process(method, radius)

    def set_all_temporal_processes(self, method=MetExtractor.DEFAULT_PROCESS_TEMPORAL,
                             period=MetExtractor.DEFAULT_PROCESS_PERIOD):
        for extractor in self._extractors:
            extractor.set_temporal_process(method, period)

    def extract_for_sites(self, sites, method=DEFAULT_METHOD, outfile_suffix=DEFAULT_OUTFILE_SUFFIX):
        print('extracting averaged daily met data for site list')

        file_base = self._filename_base.format(self._out_dir, '{}', '{}', '{}', outfile_suffix)

        datasets = []
        for extractor in self._extractors:
            settings = extractor.get_settings()
            settings['fname'] = extractor.add_outfile_suffix(outfile_suffix)
            dataset = extractor.extract_data_from_dict({'Source reference': extractor.get_source_reference()},
                                                       settings, False)
            #dataset = extractor.extract_data_from_props()
            if isinstance(dataset, list):
                # E.G. Pollen outputs a list of datasets
                datasets.extend(dataset)
                if self._verbose:
                    print('dataset output as list: {}'.format(extractor.get_source_reference()))
            else:
                # E.G. Rain, temp, wind output a single dataset
                datasets.append(dataset)
                if self._verbose:
                    print('dateset output as single: {}'.format(extractor.get_source_reference()))

        for site_key in sites:
            print('working on site: {}'.format(site_key))
            site_lat = sites[site_key][0]
            site_lon = sites[site_key][1]
            ts = Dataset({'Source reference': method,
                          'Time range': self.get_date_range(),
                          'Period': self._period,
                          'Latitude': site_lat,
                          'Longitude': site_lon,
                          'Altitude': sites[site_key][2]
                          })
            for dataset in datasets:
                ts.link(dataset)

            with open(file_base.format(site_key, date_range[0], date_range[1]), 'w') as dfile:

                dfile.write('{}, Latitude: {}, Longitude: {}\n'.format(site_key, site_lat, site_lon))

                cols_specific = ''
                for extractor in self._extractors:
                    cols_specific = ','.join((cols_specific, extractor._cols_specific))
                dfile.write(self._cols_base.format(cols_specific))

                for data in ts.values():
                    print('hello')
                    d_date = data['Linked data'][0][0]['Time']
                    d_rh = data['Linked data'][0][0]['Value']
                    d_temp = data['Linked data'][1][0]['Value']
                    (d_wspd, d_wdir) = polar(data['Linked data'][2][0]['Value'])
                    d_wdir = d_wdir * 57.29577951308232
                    if d_wdir < 0: d_wdir = d_wdir + 360.
                    d_rain = data['Linked data'][3][0]['Value']
                    #print('{}, {}, {}, {}, {}, {}\n'.format(d_date, d_rh, d_temp, d_wspd, d_wdir, d_rain))
                    dfile.write('{}, {}, {}, {}, {}, {}\n'.format(d_date, d_rh, d_temp, d_wspd, d_wdir, d_rain))

        print('finished extracting met data for pollen stations')


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

    AVAILABLE_MEASUREMENTS = ['pollen', 'rain', 'temp',  'wind', 'rel_hum', 'pressure', 'dewpoint', 'site_list']

    ### parameters

    # Measurement
    parser.add_argument("--measurements", "-m", metavar='M', type=str, nargs='+', help="measurements to be extracted. \
            Must be in (and defaults to) {}".format('[' + ", ".join([measure for measure in AVAILABLE_MEASUREMENTS]) + ']'))

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
        measurements = args.measurements
        print('Using measurements: {}'.format(measurements))
    else:
        print('No measurements provided, so using default: ', '[' + ", ".join([measure for measure in AVAILABLE_MEASUREMENTS]) + ']')
        measurements = AVAILABLE_MEASUREMENTS

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
        print('No latitude_range provided, so using default: {}'.format(str(MetExtractor.UK_LATITUDES)[1:-1].replace(',','')))
        latitude_range = MetExtractor.UK_LATITUDES

    if args.longitude_range:
        if len(args.longitude_range) >= 2:
            longitude_range = args.longitude_range[0:2]
        else:
            raise ValueError('Unable to obtain 2 values from input --longitude_range: {}'.format(str(args.longitude_range)))
        print('Using longitude range: {}'.format(str(longitude_range)[1:-1].replace(',','')))
    else:
        print('No longitude_range provided, so using default: {}'.format(str(MetExtractor.UK_LONGITUDES)[1:-1].replace(',','')))
        longitude_range = MetExtractor.UK_LONGITUDES

    if args.verbose:
        VERBOSE = max(args.verbose, 0)
        print('output verbose level: {}'.format(VERBOSE))
    else:
        print('No verbose flag provided, so using default: {}'.format(str(MetExtractor.DEFAULT_VERBOSE)))
        VERBOSE = MetExtractor.DEFAULT_VERBOSE


    # Check directory name is valid and create directory
    try:
        print('Creating directory: {}, unless it already exists.'.format(outdir_name))
        create_directory(outdir_name)
    except:
        raise ValueError('Unable to create directory: {}'.format(outdir_name))


    ### Prepare inputs and perform data extraction

    # Create necessary extractors
    if 'site_list' in measurements:
        extractors = []

        if 'pollen' in measurements:
            extractors.append(MetExtractorPollen(outdir_name))

        if 'rel_hum' in measurements:
            extractors.append(MetExtractorRelativeHumidity(outdir_name))
        if 'temp' in measurements:
            extractors.append(MetExtractorTemp(outdir_name))
        if 'wind' in measurements:
            extractors.append(MetExtractorWind(outdir_name))
        if 'rain' in measurements:
            extractors.append(MetExtractorRain(outdir_name))
        if 'pressure' in measurements:
            extractors.append(MetExtractorStationPressure(outdir_name))
        if 'dewpoint' in measurements:
            extractors.append(MetExtractorDewpoint(outdir_name))


        met_extractor_group = MetExtractorGroup(extractors, outdir_name)
        met_extractor_group.set_spatial_processes('sp_idw_mean', 50000)
        met_extractor_group.set_spatial_temperature_process('sp_altadj_idw_mean', 50000)
        met_extractor_group.set_all_temporal_processes('tp_mean', 1)
        met_extractor_group.set_all_date_ranges(date_range)
        met_extractor_group.set_all_regions(latitude_range, longitude_range)
        met_extractor_group.extract_for_sites(MetExtractorPollen.pollen_sites, 'MEDMI TS', outfile_suffix)
    else:
        if 'rain' in measurements:
            met_extractor_rain = MetExtractorRain(outdir_name)
            met_extractor_rain.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
        if 'temp' in measurements:
            met_extractor_temp = MetExtractorTemp(outdir_name)
            met_extractor_temp.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
        if 'rel_hum' in measurements:
            met_extractor_relhum = MetExtractorRelativeHumidity(outdir_name)
            met_extractor_relhum.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
        if 'pressure' in measurements:
            met_extractor_press = MetExtractorStationPressure(outdir_name)
            met_extractor_press.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
        if 'dewpoint' in measurements:
            met_extractor_dewpoint = MetExtractorDewpoint(outdir_name)
            met_extractor_dewpoint.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
        if 'wind' in measurements:
            met_extractor_wind = MetExtractorWind(outdir_name)
            met_extractor_wind.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
        if 'pollen' in measurements:
            met_extractor_pollen = MetExtractorPollen(outdir_name)
            met_extractor_pollen.extract_data(date_range, latitude_range, longitude_range, outfile_suffix)
