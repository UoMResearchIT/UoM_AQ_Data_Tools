"""
Created on Tue Sep 22 09:16:40 2020

Unified production script for MEDMI data processing. This will:
    1) load and clean our dataset
        a) find duplicated readings
            i)  filtering out the METAR data by lack of pressure reading
            ii) retaining the first value where there's no difference in presence of pressure reading
        b) removing stations identified as unwanted
            i) station 117 is on top of a mountain in the Cairngorms - RH readings are suspect,
                    and as it is unlikely to be useful comparison with participant data, we will remove it
        c) find and remove the synoptic spot readings (these are single readings per day - so we will
                     identify all such single readings and remove them, if they are synoptic spot readings or not)
        (Points (a) and (c) based on pers. comms. with Martyn Sunter, Met Office, July 2020.
         Point (b) based on data exploration by authors.)
    2)

@author: mbessdl2
"""

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

    def process(self, file_in, outfile_suffix='', date_range=None,
                exclude_site_list=DEFAULT_EXCLUDE_STATION_LIST,
                min_temperature=DEFAULT_MIN_TEMPERATURE, reference_num_stations=DEFAULT_REFERENCE_NUMBER_STATIONS,
                min_years=DEFAULT_MIN_YEARS, min_years_reference=DEFAULT_MIN_YEARS_REFERENCE,
                impute_data=PostProcessor.DEFAULT_IMPUTE_DATA, print_stats=PostProcessor.DEFAULT_PRINT_STATS,
                random_state=DEFAULT_IMPUTER_RANDOM_STATE, add_indicator=DEFAULT_IMPUTER_ADD_INDICATOR,
                initial_strategy=DEFAULT_IMPUTER_INITIAL_STRATEGY,
                max_iter=DEFAULT_IMPUTER_MAX_ITER, estimator=DEFAULT_IMPUTER_ESTIMATOR,
                output_distribution=DEFAULT_TRANSFORMER_OUTPUT_DISTRIBUTION,
                save_to_csv=PostProcessor.DEFAULT_SAVE_TO_CSV):

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
        met_extracted_data = self.load_met_data(file_in)

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
        met_data_hourly = self.combine_and_organise_mean_max(
            met_extracted_data, met_data_temp, met_data_pres, met_data_rh)

        if save_to_csv:
            # write data to file
            if self.verbose > 1: print('Writing to file: {}'.format(self.file_out.format(outfile_suffix)))
            met_data_hourly.to_csv(self.file_out, index=True, header=True, float_format='%.2f')

        return met_data_hourly


    def load_met_data(self, file_in):
        print('    load data file')
        try:
            met_data = pd.read_csv(file_in, usecols=self.get_all_columns(), skiprows=self._skip_input_rows,
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
        arguments:
            var_string:
                label for variable of interest

        other dependancies:
            met_extracted_data:
                measurement data set
            self.min_years (default 1):
                minimum number of years of data that a site must have
            self.min_years_reference (default 3.5):
                minimum number of years of data for any site that we
                are going to use as a reference site later

        returns:
            required_site_list:
                list of sites with a data count > min_years
            reference_site_list:
                list of sites with a data count > useful_num_years
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


    #%% station geographic routines

    def get_station_distances(self, site_in, useful_sites_in):

        station_location = (self.station_data.loc[site_in]['Latitude'], self.station_data.loc[site_in]['Longitude'])
        station_distances = self.calc_station_distances(stations_in=self.station_data.loc[useful_sites_in], \
                                                   stat_location=station_location)

        # sort by distance, then drop any station which is the same location as our site of interest
        station_distances = station_distances.sort_values(by='Distance',ascending=True)
        station_distances[station_distances.Distance==0]=np.nan
        station_distances = station_distances.dropna()

        return station_distances

    #%% functions for imputation of the datasets

    def transform_and_impute_data(self, df_in):
        # copy the input array, and note the columns
        df_work = df_in.copy(deep=True)
        cols = df_in.columns

        # find missing datasets to remove
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
        print('sorting temperature data')
        met_data_out_temp = self.sort_datasets(met_extracted_data, req_sites_temp, 'temperature')
        print('sorting pressure data')
        met_data_out_pressure = self.sort_datasets(met_extracted_data, req_sites_pres, 'pressure')
        print('sorting dew point temperature')
        met_data_out_dewpoint = self.sort_datasets(met_extracted_data, req_sites_dewpoint, 'dewpoint')

        return met_data_out_temp, met_data_out_pressure, met_data_out_dewpoint

    def remove_low_temperature_data(self, met_data_in):
        md_cold = met_data_in[met_data_in['temperature'] < self.min_temperature]

        if len(md_cold) > 0:
            print('     deleting this temperature and dew point temperature data')
            print(md_cold)
            met_data_in.loc[md_cold.index, ['temperature', 'dewpoint']] = np.nan
        else:
            print('    no low temperature data to remove')

        return met_data_in

    def extract_mean_max(self, ts_in, var_in_string, var_out_string):
        temp_groups = ts_in.groupby([pd.Grouper(level="date", freq='1D'), 'siteID'])
        out_data = pd.DataFrame()

        out_data['{}.max'.format(var_out_string)] = temp_groups.max()[var_in_string]
        out_data['{}.mean'.format(var_out_string)] = temp_groups.mean()[var_in_string]
        out_data['{}.flag'.format(var_out_string)] = temp_groups.mean()['{}.flag'.format(var_in_string)]

        return out_data


    def combine_and_organise_mean_max(self, met_data_in, met_data_temp, met_data_pres, met_data_rh):

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
