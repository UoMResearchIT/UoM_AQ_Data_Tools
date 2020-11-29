import pandas as pd
import numpy as np
from datetime import datetime

try:
    from sklearn.experimental import enable_iterative_imputer
    from sklearn.impute import IterativeImputer
    from sklearn.linear_model import BayesianRidge
    from sklearn import preprocessing

    import metpy.calc as mpcalc
    from metpy.units import units
except:
    pass  # print('Warning: Unable to load library: {}'.format(err))

from environmental_data_modules import PostProcessor, MetModule, DateRangeProcessor


class MetPostProcessor(PostProcessor, MetModule, DateRangeProcessor):
    """
        Class used for post-processing data that has been extracted from MEDMI server.
    """

    # Define 'absolute' constants
    BASE_FILE_OUT = '{}/Met_ppd_daily_mean_max_temp_RH_pres{}.csv'
    DATE_CALCS_FORMAT = '%Y-%m-%d %H:%M:%S'
    COLUMNS_SPECIFIC = ['temperature', 'rel_hum', 'pressure', 'dewpoint']

    # Define default constants
    DEFAULT_OUT_DIR = 'Met_processed_data'
    DEFAULT_FILE_IN = 'data_met/temp_rh_press_dewpoint_2016-2019.csv'
    DEFAULT_STATION_DATA_FILENAME = "../station_data/station_data_clean.csv"
    DEFAULT_EXCLUDE_STATION_LIST = [117]

    # Calculation defaults
    DEFAULT_MIN_TEMPERATURE = -20
    DEFAULT_REFERENCE_NUMBER_STATIONS = 5
    DEFAULT_MIN_YEARS_REFERENCE = 0.0625 #3.5
    DEFAULT_MIN_YEARS = 0.04 #1
    DEFAULT_IMPUTER_RANDOM_STATE = 0
    DEFAULT_IMPUTER_ADD_INDICATOR = True
    DEFAULT_IMPUTER_INITIAL_STRATEGY = 'mean'
    DEFAULT_IMPUTER_MAX_ITER = 300
    DEFAULT_IMPUTER_ESTIMATOR = None
    DEFAULT_TRANSFORMER_OUTPUT_DISTRIBUTION = 'normal'

    def __init__(self, out_dir=DEFAULT_OUT_DIR, station_data_filename=DEFAULT_STATION_DATA_FILENAME,
                 verbose=PostProcessor.DEFAULT_VERBOSE):
        """ Initialise instance of the MetPostProcessor class.
            Initialises the private class variables

            Args:
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of MetPostProcessor

        """
        super(MetPostProcessor, self).__init__(out_dir, verbose)
        MetModule.__init__(self)
        DateRangeProcessor.__init__(self)
        self._columns_specific = MetPostProcessor.COLUMNS_SPECIFIC
        self.station_data = station_data_filename
        self.min_temperature = MetPostProcessor.DEFAULT_MIN_TEMPERATURE
        self.min_years = MetPostProcessor.DEFAULT_MIN_YEARS
        self.min_years_reference = MetPostProcessor.DEFAULT_MIN_YEARS_REFERENCE
        self.reference_num_stations = MetPostProcessor.DEFAULT_REFERENCE_NUMBER_STATIONS
        self.impute_data = False
        self.imputer = None
        self.transformer = None

    @PostProcessor.transformer.setter
    def transformer(self, transformer):
        if transformer is None or type(transformer).__name__ == 'QuantileTransformer':
            self._transformer = transformer
        else:
            raise ValueError('Error setting transformer, incorrect object type: {}'.format(type(transformer).__name__))

    @PostProcessor.station_data.setter
    def station_data(self, filename):
        try:
            self._station_data = pd.read_csv(filename)
        except Exception as err:
            raise ValueError('Error loading station_data file {}. {}'.format(filename, err))
        else:
            try:
                self._station_data = self._station_data.set_index('Station')
            except ValueError:
                raise ValueError('Station data file has no column header: Station')

    @property
    def min_temperature(self):
        return self._min_temperature

    @min_temperature.setter
    def min_temperature(self, min_temp):
        try:
            min_years = float(min_temp)
        except ValueError:
            raise ValueError('min_temperature value ({}) must be numeric'.format(min_temp))
        # Todo: More tests to force a range of allowable temps
        self._min_temperature = min_temp

    @property
    def reference_num_stations(self):
        return self._reference_num_stations

    @reference_num_stations.setter
    def reference_num_stations(self, num_stations):
        try:
            num_stations = int(num_stations)
        except ValueError:
            raise ValueError('reference_num_stations value ({}) must be an integer'.format(num_stations))
        if num_stations < 0:
            raise ValueError('reference_num_stations value ({}) must be non-negative.')
        self._reference_num_stations = num_stations

    def process(self, in_file, outfile_suffix='', date_range=None,
                exclude_site_list=DEFAULT_EXCLUDE_STATION_LIST,
                min_temperature=DEFAULT_MIN_TEMPERATURE, reference_num_stations=DEFAULT_REFERENCE_NUMBER_STATIONS,
                min_years=DEFAULT_MIN_YEARS, min_years_reference=DEFAULT_MIN_YEARS_REFERENCE,
                impute_data=PostProcessor.DEFAULT_IMPUTE_DATA, print_stats=PostProcessor.DEFAULT_PRINT_STATS,
                random_state=DEFAULT_IMPUTER_RANDOM_STATE, add_indicator=DEFAULT_IMPUTER_ADD_INDICATOR,
                initial_strategy=DEFAULT_IMPUTER_INITIAL_STRATEGY,
                max_iter=DEFAULT_IMPUTER_MAX_ITER, estimator=DEFAULT_IMPUTER_ESTIMATOR,
                output_distribution=DEFAULT_TRANSFORMER_OUTPUT_DISTRIBUTION,
                save_to_csv=PostProcessor.DEFAULT_SAVE_TO_CSV):
        """ Post process the data extracted from MEDMI server, based on parameters

            Args:
                in_file:                (str) The file spec of the input file (required)
                date_range:             (list of 2 datetime) The date range of interest
                exclude_site_list:      (list of string/number) Site IDs to be ignored
                min_years_reference:    (float) The minimum number of years of data for any site that  we are going to
                                            use as a reference site later. (this cannot be less than min_years)
                min_years:              (float) The minimum number of years of data that a site must have
                min_temperature:        (float) The minimum temperature to be used (lower are ignored)
                reference_num_stations: (int) The number of stations to be used for imputation
                impute_data:            (boolean) Whether to attempt to impute missing data
                print_stats:            (boolean) Whether to printout the calculation statistics
                random_state:           (int) (IterativeImputer & QuantileTransformer) seed for pseudo random number generator
                output_distribution:    (str) (QuantileTransformer) Marginal distribution for the transformed data
                add_indicator:          (boolean) (IterativeImputer) if True adds a `MissingIndicator` transform to the stack
                initial_strategy:       (str) (IterativeImputer) define strategy to use for initialising missing values
                max_iter:               (int) (IterativeImputer) maximum number of imputation rounds to perform
                estimator:              (IterativeImputer) estimator method to be used
                save_to_csv:            (boolean) Whether to save the output dateframes to CSV file(s)
                outfile_suffix:         (str) The suffix to appended to the end of output file names.

            Returns:
                met_data_daily: daily dataset, for all measurements, as pandas.Dataframe
                    Required MultiIndex:
                        'time_stamp'  (datetime object): date (only) (e.g. 2017-06-01)
                        'sensor_name'          (string): ID string for site (e.g. '3 [WEATHER]')
                    Required columns:
                        'Temperature.max'       (float): daily maximum value
                        'Temperature.mean'      (float): daily mean value
                        'Temperature.flag'      (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
                        'RelativeHumidity.max'  (float): daily maximum value
                        'RelativeHumidity.mean' (float): daily mean value
                        'RelativeHumidity.flag' (float): flag to indicate fraction of imputed data
                                                            (1 = fully imputed, 0 = no imputed values were used)
                        'Pressure.max'          (float): daily maximum value
                        'Pressure.mean'         (float): daily mean value
                        'Pressure.flag'         (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
        """

        if date_range is not None:
            self.date_range = [datetime.strptime(date_range[0], MetPostProcessor.INPUT_DATE_FORMAT),
                               datetime.strptime(date_range[1], MetPostProcessor.INPUT_DATE_FORMAT)]
        else:
            self.date_range = [self.get_available_start(), self.get_available_end()]

        self.min_temperature = min_temperature
        self.reference_num_stations = reference_num_stations
        self.min_years = min_years
        self.min_years_reference = min_years_reference
        self.file_out = MetPostProcessor.BASE_FILE_OUT.format(self.out_dir, outfile_suffix)
        self._print_stats = print_stats
        if impute_data:
            self.imputer = IterativeImputer(random_state=random_state, add_indicator=add_indicator,
                                       initial_strategy=initial_strategy, max_iter=max_iter, verbose=self.verbose,
                                       estimator=estimator)
        self.impute_data = impute_data

        # set the power transform options
        self.transformer = preprocessing.QuantileTransformer(output_distribution=output_distribution,
                                                             random_state=random_state)

        print('checking validity of and loading met data file')
        met_extracted_data = self.load_met_data(in_file)

        if self.verbose > 1: print('Metadata before dropping/filtering: \n {}'.format(met_extracted_data))

        print('dropping duplicate values and unwanted stations')
        met_extracted_data = self.find_and_drop_duplicates_and_unwanted_stations(met_extracted_data,
                                                                                 exclude_site_list)

        print('dropping single daily measurement stations')
        met_extracted_data = self.drop_single_daily_measurement_stations(met_extracted_data)

        print('filtering to remove unrealistically low temperatures')
        met_extracted_data = self.remove_low_temperature_data(met_extracted_data)

        print('filter for minimum data lengths, and reduce dataset to only stations of interest')
        met_extracted_data, reference_sites, req_sites_temp, req_sites_pres, req_sites_dewpoint = \
            self.list_required_and_reference_sites(met_extracted_data)

        if self.verbose > 1: print('Metadata after dropping/filtering: {}'.format(met_extracted_data))

        if len(met_extracted_data.index) == 0:
            print('Exiting post-processing: Metadata is empty after initial filtering processes')
            return

        if self.impute_data:
            print('imputation of data, returning hourly data')
            met_data_temp, met_data_pres, met_data_dewpoint = self.organise_data_imputation(
                met_extracted_data, reference_sites, req_sites_temp, req_sites_pres, req_sites_dewpoint)
        else:
            print('sorting data (no imputation), returning hourly data')
            met_data_temp, met_data_pres, met_data_dewpoint = self.organise_data(
                met_extracted_data, req_sites_temp, req_sites_pres, req_sites_dewpoint)

        print('calculation of relative humidity from temperature and dew point temperature')
        met_data_rh = self.rh_calculations(met_extracted_data, met_data_temp, met_data_dewpoint)

        # calculate the daily max and mean for each station
        met_data_daily = self.combine_and_organise_mean_max(
            met_extracted_data, met_data_temp, met_data_pres, met_data_rh)

        if save_to_csv:
            # write data to file
            if self.verbose > 1: print('Writing to file: {}'.format(self.file_out.format(outfile_suffix)))
            met_data_daily.to_csv(self.file_out, index=True, header=True, float_format='%.2f')

        return met_data_daily


    def load_met_data(self, file_in):
        '''
        Loading the meteorological dataset.
        
        Args:
            file_in (Path object): path for the file that is to be read in
        
        Uses:
            self.get_all_column_headers
            self._skip_input_rows
        
        Returns:
            met_data: Pandas dataframe containing meteorological dataset. 
                      This will include corrected date strings.
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
        '''
    
        print('    load data file')
        try:
            met_data = pd.read_csv(file_in, usecols=self.get_all_column_headers(), skiprows=self._skip_input_rows,
                                   engine='python', sep=',\s*', na_values='None')
        except Exception as err:
            raise ValueError('Error reading specified file: {}. {}'.format(file_in, err))
        if len(met_data.index) < 1:
            raise ValueError('Input file ({}) is empty'.format(file_in))
        print('    correct date string')
        met_data['date'] = met_data['date'].apply(self.parse_calcs_date)

        return met_data


    #%% function for writing out some information about the data count stats

    def print_data_count_stats(self, dc_in):
        '''
        This prints the daily data count statistics. It will output the total number
        of days with data (this should be equal to full time period requested for the
        dataset), then print out the count of days containing a set number of data 
        values (between 0 and 48) for each of the specified variables.
        
        Args:
            dc_in: Pandas groupby dataframe containing met data, 
                    grouped by site and date (at 1 day frequency)
                Required index:
                    'siteID'      (string) ID string for site
                    'date'        (datetime object) date of measurement set
                Required columns:
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature

        Uses:
            self.COLUMNS_SPECIFIC: list of strings for selecting columns to output 
                                   (and the order for these to be output!)
        '''

        print('total temperature daily data count is: {}'.format(dc_in.count().values[0]))
        print('# data points per day, total daily data point counts')
        print('      temperature, rel hum, pressure, dew point temp')

        for dpoint in range(0,49):
            dcounts = [dpoint]
            for col_name in self.COLUMNS_SPECIFIC:
                dcounts.append(dc_in[dc_in[col_name]==dpoint].count().values[0])
            print('{0:4},{1:7},{2:7},{3:7},{4:7}'.format(*dcounts))

    #%% functions for finding and removing unwanted data

    def find_and_drop_duplicates_and_unwanted_stations(self, met_data_in, exlude_stations):
        '''
        This cleans the meteorological dataset, to remove duplicated measurements, and
        to remove unwanted stations from the dataset.
        
        Args:
            met_data_in: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature

            exlude_stations: (list, strings) station numbers to exclude from final dataset
        
        Returns:
            met_data_in: cleaned pandas dataframe, with same data structure as 'met_data_in'
            
        '''
        # Pull out all duplicated values
        met_duplicates = met_data_in[met_data_in.duplicated(subset=['date', 'siteID'], keep=False)]

        # Split these into those with, and without, pressure data (SYNOP will have pressure data)
        #   We will keep (most of) the readings with pressure data, and will drop (most of) the
        #   readings without pressure data.
        met_dup_with_pres = met_duplicates[met_duplicates['pressure'].isna() == False]
        met_dup_no_pres = met_duplicates[met_duplicates['pressure'].isna() == True]

        # Get the 2nd data points which are duplicated *and* where both have pressure data (these are only a few points)
        #   These will be the few readings with pressure data that we drop, so we set keep to 'first' so that we
        #   get the indexes for the 'last' values.
        met_dup_with_pres_dups = met_dup_with_pres[met_dup_with_pres.duplicated(subset=['date', 'siteID'], keep='first')]

        # Find the duplicates where there is no pressure data for either, and save the first value (<10,000 points)
        #   In this case we also want to preserve the first data points, but as we are going to drop the
        #   indexes we extract from our list of indexes to throw away we want to set keep to 'last' in this instance.
        met_dup_no_pres_dups = met_dup_no_pres[met_dup_no_pres.duplicated(subset=['date', 'siteID'], keep='last')]
        met_dup_no_pres_reduced = met_dup_no_pres.drop(index=met_dup_no_pres_dups.index)


        # Build an index of the duplicated datapoints to drop - starting with the data with no pressure readings
        #   that will be dropped, then appending the indexes of the few data points with pressure data that we don't want.
        indexes_to_drop = met_dup_with_pres_dups.index
        indexes_to_drop = indexes_to_drop.append(met_dup_no_pres_reduced.index)


        # Append to this list the indexes of all station data that we are dropping completely
        for station in exlude_stations:
            station_drop_indexes = met_data_in[met_data_in['siteID'] == station].index
            indexes_to_drop = indexes_to_drop.append(station_drop_indexes)

        # Finally drop all the data that is unwanted
        return met_data_in.drop(index=indexes_to_drop)

    def drop_single_daily_measurement_stations(self, met_data_in):
        '''
        Searches the dataset to find stations which *only* have single measurements each
        daily, and removes these stations. Stations which have a mix of single daily 
        measurements, and days with multiple measurements, are kept.
        
        Args:
            met_data_in: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            
        Uses:
            self._print_stats (logical)
        
        Returns:
            met_data_reduced: cleaned pandas dataframe, with same data structure as 'met_data_in'
        
        '''
        # group the data by date, and count the readings per day
        tempgroups = met_data_in.groupby(['siteID', pd.Grouper(key='date', freq='1D')])
        data_counts = tempgroups.count()

        # some diagnostic output, if required
        if self._print_stats:
            self.print_data_count_stats(data_counts)

        # find the stations and days with single temperature measurements for that day
        temperature1_index = data_counts[data_counts['temperature']==1].index
        # find the stations and days with more than a single temperature measurement
        temperature24_index = data_counts[data_counts['temperature']>1].index

        # get the station ID's for both of these indexes
        station1_list  = temperature1_index.get_level_values(0).unique().to_list()
        station24_list = temperature24_index.get_level_values(0).unique().to_list()

        # Identify the stations in list 1 but not list 24 - these will be the stations with only single daily readings
        #   Note: A few stations which seem to be mixed, with long(ish) periods of single daily readings,
        #         mixed with longer periods of multiple daily readings will be missed by this rough filter.
        #         We will be filtering later by total hourly data counts, which will catch any examples of
        #         this where the data sets are too sparse to be used.
        station_single_list = [x for x in station1_list if x not in station24_list]

        if self._print_stats:
            print('single daily measurement stations to drop')
            print(station_single_list)

        # drop all the stations with only single daily readings
        indexes_to_drop = pd.Int64Index(data=[],dtype='int64')
        for station in station_single_list:
            station_drop_indexes = met_data_in[met_data_in['siteID']==station].index
            indexes_to_drop = indexes_to_drop.append(station_drop_indexes)

        # drop all the data that is unwanted
        met_data_reduced = met_data_in.drop(index=indexes_to_drop)

        # print diag output for new dataset
        if self._print_stats:
            print('new stats for the reduced data:')
            temp_groups = met_data_reduced.groupby(['siteID', pd.Grouper(key='date', freq='1D')])
            data_counts = temp_groups.count()
            self.print_data_count_stats(data_counts)

        return met_data_reduced


    #%% function for getting two lists of stations, one for required site, one for reference sites

    def station_listing(self, met_extracted_data, var_string):
        #Todo Doug: can this be merged with station_listing() in aurn_post_processor (and put in to post_processor)
        '''
        This function lists the stations which fulfil two requirements:
            1) those that meet minimum data count for inclusion in the dataset
            2) those that meet minimum data count for use as a reference station
        
        Args:
            var_string (string): label for variable of interest, used to 
                                 select column containing data of interest
            met_extracted_data: hourly met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature

        Uses:
            self.min_years (float, default 1):
                minimum number of years of data that a site must have to be included 
                in the final dataset.
            self.min_years_reference (float, default 3.5):
                minimum number of years of data for any site that is going to be used 
                as a reference site.

        Returns:
            required_site_list (list, string or int):
                                 sites with a data count > min_years
            reference_site_list (list, string or int):
                                 sites with a data count > useful_num_years
        '''

        site_list_interior = met_extracted_data['siteID'].unique()

        required_site_list = []
        reference_site_list = []

        for site in site_list_interior:
            metsite = met_extracted_data[met_extracted_data['siteID'] == site]
            try:
                measurement_num = len(metsite[metsite[var_string].notna()])
            except:
                measurement_num = 0
            if measurement_num > self.min_years*365*24:
                required_site_list.append(site)
            if measurement_num > self.min_years_reference*365*24:
                reference_site_list.append(site)

        return required_site_list, reference_site_list


    def list_required_and_reference_sites(self, met_data_in):
        '''
        This function creates the lists of required sites, and reference sites, for the 
        final dataset.
        
        Args:
            met_data_in: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            
        Returns:
            met_data_filtered: pandas dataframe, as above, containing hourly dataset for only 
                               the required station datasets
            reference_sites: (list, string or int) the siteID's for our reference sites 
                             (composite list made up from the measurement specific lists below)
            req_sites_temp: (list, string or int) required sites for temperature data
            req_sites_pres: (list, string or int) required sites for pressure data
            req_sites_dewpoint: (list, string or int) required sites for dewpoint temperature data
        
        '''
        print('    get the lists of required and reference stations for each measurement variable')
        req_sites_temp, reference_sites_temp = self.station_listing(met_data_in, 'temperature')
        req_sites_pres, reference_sites_pres = self.station_listing(met_data_in, 'pressure')
        req_sites_dewpoint, reference_sites_dewpoint = self.station_listing(met_data_in, 'dewpoint')

        # find a unified list of useful sites for all our measurements
        reference_sites = [x for x in reference_sites_dewpoint if x in reference_sites_pres]

        # get a list of all sites which are required for at least one measurement set
        required_sites = list(dict.fromkeys(req_sites_temp + req_sites_pres + req_sites_dewpoint))
        print('there are {} required sites, and {} reference sites'.format(len(required_sites), len(reference_sites)))
        met_data_filtered = met_data_in[met_data_in['siteID'].isin(required_sites)]

        return met_data_filtered, reference_sites, req_sites_temp, req_sites_pres, req_sites_dewpoint

    #%% functions for calculating RH from temperature and dew point temperature data

    def rh_calculations(self, met_data_in, met_data_temp, met_data_dewpoint):
        '''
        This function calculates the relative humidity data, using the metpy
        relative_humidity_from_dewpoint function.
        
        Args:
            met_data_in: met data as a pandas.DataFrame, used for original RH dataset
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            met_data_temp: temperature data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'temperature'    (float): temperature
                    'temperature.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
            met_data_dewpoint: dewpoint temperature data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'dewpoint'    (float): dewpoint temperature
                    'dewpoint.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
            
        Uses:
            self.verbose (int)
            self._print_stats (logical)
            
        Returns:
            met_data_out: relative humidity data as pandas.Dataframe 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'rel_hum'    (float): relative humidity
                    'rel_hum.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed),
                                         data is marked as imputed if at least one of temperature or
                                         relative humidity has been imputed.
        '''
        # merge the two input datasets, dropping indexes which are not in both
        met_data_all = met_data_temp.merge(met_data_dewpoint, how='inner', left_index=True, right_index=True)
        if self.verbose > 1: print('Met data in: {}'.format(met_data_all))

        # create output data frame
        met_data_out = pd.DataFrame(index=met_data_all.index)
        met_data_out['rel_hum'] = mpcalc.relative_humidity_from_dewpoint(
            met_data_all['temperature'].values * units.degC,
            met_data_all['dewpoint'].values * units.degC) * 100.0

        met_data_out['rel_hum.flag'] = met_data_all[['temperature.flag','dewpoint.flag']].values.max(1)

        # plot the distribution of the calculated RH difference from measured RH
        if self._print_stats:
            met_data_internal = met_data_in.set_index(['date', 'siteID'])

            met_data_internal['rh2'] = met_data_out['rel_hum']

            met_data_internal['rhdiff'] = met_data_internal['rh2'] - met_data_internal['rel_hum']

            ax = met_data_internal['rhdiff'].hist(bins=100)
            ax.semilogy()
            ax.set_xlabel('RH_new - RH')

        return met_data_out


    #%% functions for imputation of the datasets

    def transform_and_impute_data(self, df_in):
        '''
        
        Function for organising the transformation and imputation of the datasets.
        The input dataset is processed to remove missing variables (e.g. where a
        station has temperature data, but no pressure data, the pressure column will be
        removed before imputation) before transformation & imputation. Afterwards the 
        returned dataset has the missing variables readded, so that it has the same
        shape as the original dataset.
        
        Args:
            df_in: pandas dataframe, containing hourly datasets for imputation
                Required Index:
                    'date'        (datetime object) date/time of measurement
                Required columns:
                    '[var_string]'  (float): met data for the site of interest
                Optional columns (repeated [self.reference_num_stations] times):
                    '[station_list_string]' (float): met data for selected reference sites
                    
        Uses:
            self.transformer
            self.imputer
        
        Returns:
            df_out: pandas dataframe, containing imputed dataset, same shape as 
                    input dataframe. Missing data for all stations will be imputed
                    and returned, but only the data for '[var_string]' is retained
                    in the end. 
        '''
        # copy the input array, and note the columns
        df_work = df_in.copy(deep=True)
        cols = df_in.columns

        # find missing datasets to remove (if any)
        # also we note the columns that will be saved, and their order, for transferring data back!
        col_remove = []
        col_save = []
        for col in cols:
            if all(df_work[col].isna()):
                col_remove.append(col)
            else:
                col_save.append(col)
        df_work = df_work.drop(columns=col_remove)

        if self.verbose > 1: print('df_work input to quantile transformer: \n {}'.format(df_work))

        # power transformer fitting and transforming
        self.transformer.fit(df_work.dropna())
        np_out = self.transformer.transform(df_work)

        # impute the missing values in this new dataframe
        self.imputer.fit(np_out)
        imp_out = self.imputer.transform(np_out)

        # apply the inverse transformation for our datasets (leaving out the indicator flags)
        np_inv = self.transformer.inverse_transform(imp_out[:, :np_out.shape[1]])

        # copy the transformed values to a new dataframe
        df_out = df_in.copy(deep=True)
        for pos,col in enumerate(col_save):
            pos_out = list(cols).index(col)
            df_out.iloc[:,pos_out] = np_inv[:,pos]

        return df_out

    def get_full_datasets(self, met_extracted_data, req_sites_list, useful_sites_list, station_list_string, var_string):
        '''
        Function for the imputation and building of the full hourly datasets.
        
        Args:
            met_extracted_data: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            
            req_sites_list (list, string or int): sites to retain for data
            useful_sites_list (list, string or int): sites to be used as reference sites
            station_list_string (list, string or int): strings used to identify 
                                 reference stations in the working dataset (not site specific)
            var_string     (string): name of the variable to sort dataset for
            
        Uses:
            self.start
            self.end
            self.reference_num_stations
        
        Returns:
            full_data_out: selected meteorological data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    '[var_string]'    (float): met variable data
                    '[var_string].flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
        '''

        date_index = pd.date_range(start=self.start, end=self.end,
                                   freq='1H', name='date')

        # add the Date index
        indexed_orig_data = met_extracted_data.set_index('date')

        # define initial column for dataframe
        dataframe_columns = {var_string: np.nan}

        # empty dataframe for storing data
        full_data_out = pd.DataFrame()

        for site in req_sites_list:

            print('working on site {}'.format(site))

            work_data = indexed_orig_data[indexed_orig_data.siteID==site]

            ts = pd.DataFrame(dataframe_columns, index=date_index)

            ts[var_string] = work_data[var_string]

            if(len(work_data)<len(ts)):

                print('  site is missing {} data points, filling these in'.format(len(ts)-len(work_data)))

                station_distances = self.get_station_distances(site, useful_sites_list)
                # get data for the [reference_station_number] closest stations:
                for ii in range(0, self.reference_num_stations):
                    ts[station_list_string[ii]] = \
                        indexed_orig_data[indexed_orig_data.siteID==station_distances.index[ii]][var_string]

                # run the imputation process
                imputed_hourly_dataframe = self.transform_and_impute_data(ts)

                # drop the extra columns of station data
                ts = ts.drop(columns=station_list_string)

                # add a flag to denote the values which have been imputed
                #   then replace the original data with the imputed data
                ts['{}.flag'.format(var_string)] = ts[var_string].isna() * 1
                ts[var_string] = imputed_hourly_dataframe[var_string]
            else:
                # if we didn't impute anything, add zero value flags
                ts['{}.flag'.format(var_string)] = 0

            # add the site ID, and reindex
            ts['siteID'] = site
            ts = ts.reset_index().set_index(['date','siteID'])

            # copy data to new array
            full_data_out = full_data_out.append(ts)

        return full_data_out


    def organise_data_imputation(self, met_extracted_data, reference_sites, req_sites_temp, req_sites_pres,
                                 req_sites_dewpoint):
        '''
        Function for organising the imputation of the datasets. This runs the 
        'get_full_datasets' function for each of the variables of interest.
        
        Args:
            met_extracted_data: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            
            reference_sites (list, string or int): sites to use for reference when imputing datasets
            req_sites_temp (list, string or int): sites to retain for temperature data
            req_sites_pres (list, string or int): sites to retain for pressure data
            req_sites_dewpoint (list, string or int): sites to retain for dewpoint temperature data
        
        Uses:
            self.reference_num_stations
        
        Returns:
            met_data_out_temp: pandas dataframe, containing temperature data and flag to 
                               indicate imputed data
            met_data_out_pressure: pandas dataframe, containing pressure data and flag to 
                                   indicate imputed data
            met_data_out_dewpoint: pandas dataframe, containing dewpoint temperature data 
                                   and flag to indicate imputed data
        '''


        station_list = ['station{}'.format(x+1) for x in range(0, self.reference_num_stations)]

        print('imputing temperature data')
        met_data_out_temp = self.get_full_datasets(
            met_extracted_data, req_sites_temp, reference_sites, station_list, 'temperature')

        print('imputing pressure data')
        met_data_out_pressure = self.get_full_datasets(
            met_extracted_data, req_sites_pres, reference_sites, station_list, 'pressure')

        print('imputing dew point temperature')
        met_data_out_dewpoint = self.get_full_datasets(
            met_extracted_data, req_sites_dewpoint, reference_sites, station_list, 'dewpoint')

        return met_data_out_temp, met_data_out_pressure, met_data_out_dewpoint

    def sort_datasets(self, met_extracted_data, req_sites_list, var_string):
        '''
        Function for creating the hourly dataset when we are not imputing any data.
        
        Args:
            met_extracted_data: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            
            req_sites_list (list, string or int): sites to retain for data
            var_string     (string): name of the variable to sort dataset for
        
        Uses:
            self.start
            self.end
        
        Returns:
            full_data_out: selected meteorological data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    '[var_string]'    (float): met variable data
                    '[var_string].flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
        '''
        # AG: Trim date index to be only those available in dataset: Or memory overloads and kills process.
        # Todo Doug: check OK.
        start_date = max(met_extracted_data['date'].min(), self.start)
        end_date = min(met_extracted_data['date'].max(), self.end)

        print('Start date: {}'.format(datetime.strftime(start_date, MetPostProcessor.INPUT_DATE_FORMAT)))
        print('End date: {}'.format(datetime.strftime(end_date, MetPostProcessor.INPUT_DATE_FORMAT)))
        if self.verbose > 1: print('Using date range in sort_datasets: {} to {}'.format(
            datetime.strftime(start_date, MetPostProcessor.INPUT_DATE_FORMAT),
            datetime.strftime(end_date, MetPostProcessor.INPUT_DATE_FORMAT)
        ))

        date_index = pd.date_range(start=start_date, end=end_date, freq='1H', name='date')

        # add the Date index
        indexed_orig_data = met_extracted_data.set_index('date')

        # define initial column for dataframe
        dataframe_columns = {var_string: np.nan}

        # empty dataframe for storing data
        full_data_out = pd.DataFrame()

        for site in req_sites_list:

            print('extract site {}'.format(site))

            work_data = indexed_orig_data[indexed_orig_data.siteID==site]

            ts = pd.DataFrame(dataframe_columns, index=date_index)

            ts[var_string] = work_data[var_string]

            # as we didn't impute anything, add zero value flags
            ts['{}.flag'.format(var_string)] = 0

            # add the site ID, and reindex
            ts['siteID'] = site
            ts = ts.reset_index().set_index(['date','siteID'])

            # copy data to new array
            full_data_out = full_data_out.append(ts)

        return full_data_out

    def organise_data(self, met_extracted_data, req_sites_temp, req_sites_pres, req_sites_dewpoint):
        '''
        Function for organising the creation of the datasets when no imputation is involved. 
        This runs the 'sort_datasets' function for each of the variables of interest.
        
        Args:
            met_extracted_data: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            
            req_sites_temp (list, string or int): sites to retain for temperature data
            req_sites_pres (list, string or int): sites to retain for pressure data
            req_sites_dewpoint (list, string or int): sites to retain for dewpoint temperature data
                
        Returns:
            met_data_out_temp: temperature data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'temperature'    (float): temperature
                    'temperature.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
            
            met_data_out_pressure: pressure data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'pressure'    (float): pressure
                    'pressure.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
            
            met_data_out_dewpoint: dewpoint temperature data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'dewpoint'    (float): dewpoint temperature
                    'dewpoint.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
        '''
        print('sorting temperature data')
        met_data_out_temp = self.sort_datasets(met_extracted_data, req_sites_temp, 'temperature')
        print('sorting pressure data')
        met_data_out_pressure = self.sort_datasets(met_extracted_data, req_sites_pres, 'pressure')
        print('sorting dew point temperature')
        met_data_out_dewpoint = self.sort_datasets(met_extracted_data, req_sites_dewpoint, 'dewpoint')

        return met_data_out_temp, met_data_out_pressure, met_data_out_dewpoint

    def remove_low_temperature_data(self, met_data_in):
        '''
        Function for removing datapoints with very low temperature data. Removes both 
        temperature and dew point temperature data where the temperature data is below 
        the minimum value.
        
        Args:
            met_data_in: met data as a pandas.DataFrame
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
        
        Uses:
            self.min_temperature (float): minimum temperature to use
        
        Returns:
            met_data_in: as above, without the removed data
        '''
        
        
        md_cold = met_data_in[met_data_in['temperature'] < self.min_temperature]

        if len(md_cold) > 0:
            print('     deleting this temperature and dew point temperature data')
            print(md_cold)
            met_data_in.loc[md_cold.index, ['temperature', 'dewpoint']] = np.nan
        else:
            print('    no low temperature data to remove')

        return met_data_in

    def extract_mean_max(self, ts_in, var_in_string, var_out_string):
        '''
        Function for calculating the daily mean and maximum values, as well as the
        fractional flag variable (indicating the number of imputed datapoints in that day).
        
        Designed for reading in hourly data. Could, in theory, be used for any frequency
        of measurement, and still return mean and max, but will require data with a regular
        sampling period.
        
        Args:
            ts_in: pandas dataframe, containing hourly data for the given variable
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    '[var_in_string]'    (float): met variable
                    '[var_in_string].flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
            
            var_in_string  (string): input variable name (e.g. 'temperature')
            var_out_string (string): output variable name (e.g. 'Temperature')
            
        Returns:
            out_data: pandas dataframe, containing daily mean and max data for the variable,
                      as well as the flag indicating fraction of data which is imputed (0-1).
                Required MultiIndex:
                    'time_stamp'   (datetime object): date (only) (e.g. 2017-06-01)
                    'sensor_name'           (string): ID string for site (e.g. '3 [WEATHER]')
                Required columns:
                    '[var_out_string].max'     (float): daily maximum value
                    '[var_out_string].mean'    (float): daily mean value
                    '[var_out_string].flag'    (float): flag to indicate fraction of imputed data 
                                                   (1 = fully imputed, 0 = no imputed values were used)
        '''
        
        
        temp_groups = ts_in.groupby([pd.Grouper(level="date", freq='1D'), 'siteID'])
        out_data = pd.DataFrame()

        out_data['{}.max'.format(var_out_string)] = temp_groups.max()[var_in_string]
        out_data['{}.mean'.format(var_out_string)] = temp_groups.mean()[var_in_string]
        out_data['{}.flag'.format(var_out_string)] = temp_groups.mean()['{}.flag'.format(var_in_string)]

        return out_data


    def combine_and_organise_mean_max(self, met_data_in, met_data_temp, met_data_pres, met_data_rh):
        '''
        Function for organising the calculation of daily mean, max, and flag data.
        
        Also corrects the variable names, and the station ID's, to match what is expected
        for the final dataset.
        
        Args:
            met_data_in: not used?
            met_data_temp: hourly temperature data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'temperature'    (float): temperature
                    'temperature.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
            
            met_data_pressure: hourly pressure data as pandas.Dataframe, 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'pressure'    (float): pressure
                    'pressure.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed)
            met_data_rh: hourly relative humidity data as pandas.Dataframe 
                Required MultiIndex:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                Required columns:
                    'rel_hum'    (float): relative humidity
                    'rel_hum.flag' (int): flag to indicate imputed data (1 = imputed, 0 = not imputed),

        
        Returns:
            combined_data: daily dataset, for all measurements, as pandas.Dataframe
                Required MultiIndex:
                    'time_stamp'   (datetime object): date (only) (e.g. 2017-06-01)
                    'sensor_name'           (string): ID string for site (e.g. '3 [WEATHER]')
                Required columns:
                    'Temperature.max'     (float): daily maximum value
                    'Temperature.mean'    (float): daily mean value
                    'Temperature.flag'    (float): flag to indicate fraction of imputed data 
                                                   (1 = fully imputed, 0 = no imputed values were used)
                    'RelativeHumidity.max'     (float): daily maximum value
                    'RelativeHumidity.mean'    (float): daily mean value
                    'RelativeHumidity.flag'    (float): flag to indicate fraction of imputed data 
                                                        (1 = fully imputed, 0 = no imputed values were used)
                    'Pressure.max'     (float): daily maximum value
                    'Pressure.mean'    (float): daily mean value
                    'Pressure.flag'    (float): flag to indicate fraction of imputed data 
                                                (1 = fully imputed, 0 = no imputed values were used)
                                        
            
        
        '''

        met_groups_rh = self.extract_mean_max(met_data_rh, 'rel_hum', 'RelativeHumidity')
        met_groups_temp = self.extract_mean_max(met_data_temp, 'temperature', 'Temperature')
        met_groups_pres = self.extract_mean_max(met_data_pres, 'pressure', 'Pressure')

        combined_data = met_groups_temp.merge(met_groups_rh, how='outer', left_index=True, right_index=True)
        combined_data = combined_data.merge(met_groups_pres, how='outer', left_index=True, right_index=True)

        combined_data.sort_index(level=1,inplace=True)
        combined_data.index = combined_data.index.set_levels(
            ['{} [WEATHER]'.format(x) for x in combined_data.index.levels[1]], level=1)
        combined_data.index.rename(['time_stamp', 'sensor_name'], inplace=True)

        return combined_data
