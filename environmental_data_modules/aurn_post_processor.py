"""

Script for automated downloading of AURN data for a given date range.


Summary of how many sites we will be imputing for data for each species,
and how many useful site there are for each.

NOXasNO2 has 157 req sites and 118 useful sites
O3 has 75 req sites and 66 useful sites
NO2 has 157 req sites and 118 useful sites
PM2.5 has 81 req sites and 55 useful sites
PM10 has 77 req sites and 51 useful sites
SO2 has 27 req sites and 19 useful sites

We will only impute data for sites with > 1 year of data, and only use
sites with >3.5 years of data for the imputation inputs


Full dataset information on NaN's and negative/zero numbers:
O3 has 2507061 positive values
O3 has 3139956 NaNs
O3 has 735 negative or zero values that will be replaced with NaNs
PM10 has 2287338 positive values
PM10 has 3352010 NaNs
PM10 has 8404 negative or zero values that will be replaced with NaNs
PM2.5 has 2336213 positive values
PM2.5 has 3280461 NaNs
PM2.5 has 31078 negative or zero values that will be replaced with NaNs
NO2 has 4936512 positive values
NO2 has 707558 NaNs
NO2 has 3682 negative or zero values that will be replaced with NaNs
NOXasNO2 has 4939775 positive values
NOXasNO2 has 707425 NaNs
NOXasNO2 has 552 negative or zero values that will be replaced with NaNs
SO2 has 837779 positive values
SO2 has 4806980 NaNs
SO2 has 2993 negative or zero values that will be replaced with NaNs
(note, the NaN count will include sites that we will not be imputing with data)

"""

try:
    import datetime
    import pandas as pd
    import numpy as np
    from pathlib import Path

    from sklearn.experimental import enable_iterative_imputer
    from sklearn.impute import IterativeImputer
    from sklearn.linear_model import BayesianRidge
    from sklearn import preprocessing
except:
    pass

from environmental_data_modules import PostProcessor, AurnModule, DateRangeProcessor


class AurnPostProcessor(PostProcessor, AurnModule, DateRangeProcessor):
    # Define 'absolute' constants
    BASE_FILE_OUT = '{}/aurn_processed_daily_{}.csv'

    # Define default constants
    DEFAULT_OUT_DIR = 'Aurn_processed_data'
    DEFAULT_EMEP_FILENAME = None

    # Calculation defaults
    DEFAULT_MIN_YEARS_REFERENCE = 1
    DEFAULT_MIN_YEARS = 1
    DEFAULT_IMPUTER_RANDOM_STATE = 0
    DEFAULT_IMPUTER_ADD_INDICATOR = False
    DEFAULT_IMPUTER_INITIAL_STRATEGY = 'mean'
    DEFAULT_IMPUTER_MAX_ITER = 100
    try:
        DEFAULT_IMPUTER_ESTIMATOR = BayesianRidge()
    except:
        DEFAULT_IMPUTER_ESTIMATOR = None
    DEFAULT_TRANSFORMER_METHOD = 'box-cox'
    DEFAULT_TRANSFORMER_STANDARDIZE = False

    def __init__(self, metadata_filename=AurnModule.DEFAULT_METADATA_FILE, metadata_url=AurnModule.DEFAULT_METADATA_URL,
                 out_dir=DEFAULT_OUT_DIR, verbose=PostProcessor.DEFAULT_VERBOSE):
        super(AurnPostProcessor, self).__init__(out_dir, verbose)
        AurnModule.__init__(self, metadata_filename=metadata_filename, metadata_url=metadata_url)
        DateRangeProcessor.__init__(self)

        self._emep_data = None
        self.min_years_reference = AurnPostProcessor.DEFAULT_MIN_YEARS_REFERENCE
        self.min_years = AurnPostProcessor.DEFAULT_MIN_YEARS
        self.impute_data = False
        self.imputer = None
        self.transformer = None

    @PostProcessor.transformer.setter
    def transformer(self, transformer):
        if transformer is None or type(transformer).__name__ == 'PowerTransformer':
            self._transformer = transformer
        else:
            raise ValueError('Error setting transformer, incorrect object type: {}'.format(type(transformer).__name__))

    @PostProcessor.station_data.setter
    def station_data(self, raw_data):
        if self.verbose > 0:
            print('Loading stations data metadata')
        try:
            station_data = raw_data.drop_duplicates()
            station_data = station_data.rename(
                columns={"site_id": "SiteID", "latitude": "Latitude", "longitude": "Longitude"})
            station_data = station_data.set_index('SiteID')
        except Exception as err:
            raise ValueError('Unable to get correct site data from Metadata input file. Check metadata file content.')

        self._station_data = station_data

    def process(self, in_file, date_range=None,
                site_list=AurnModule.DEFAULT_SITE_LIST,
                emep_filename=DEFAULT_EMEP_FILENAME,
                min_years_reference=DEFAULT_MIN_YEARS_REFERENCE, min_years=DEFAULT_MIN_YEARS,
                impute_data=PostProcessor.DEFAULT_IMPUTE_DATA,
                random_state=DEFAULT_IMPUTER_RANDOM_STATE, add_indicator=DEFAULT_IMPUTER_ADD_INDICATOR,
                initial_strategy=DEFAULT_IMPUTER_INITIAL_STRATEGY,
                max_iter=DEFAULT_IMPUTER_MAX_ITER, estimator=DEFAULT_IMPUTER_ESTIMATOR,
                transformer_method=DEFAULT_TRANSFORMER_METHOD, transformer_standardize=DEFAULT_TRANSFORMER_STANDARDIZE,
                save_to_csv=PostProcessor.DEFAULT_SAVE_TO_CSV,
                outfile_suffix=''):

        # Process inputs
        if date_range is not None:
            self.date_range = [datetime.strptime(date_range[0], DateRangeProcessor.INPUT_DATE_FORMAT),
                               datetime.strptime(date_range[1], DateRangeProcessor.INPUT_DATE_FORMAT)]
        else:
            self.date_range = [self.get_available_start(), self.get_available_end()]

        self.file_out = AurnPostProcessor.BASE_FILE_OUT.format(self.out_dir, outfile_suffix)
        self._emep_data = self.load_emep_data(emep_filename)
        self.min_years = min_years
        self.min_years_reference = min_years_reference
        self.site_list = site_list
        self.station_data = self.metadata['AURN_metadata'][['site_id', 'latitude', 'longitude', 'site_name']]
        if self.verbose > 1: print('Station data: \n {}'.format(self.station_data))

        # Read in hourly dataframe file
        try:
            hourly_dataframe = pd.read_csv(in_file,
                                           sep=',',
                                           usecols=[AurnModule.INDEX_EXTRACTED].append(AurnModule.EXTRACTED_FILE_COLS),
                                           index_col=AurnModule.INDEX_EXTRACTED,
                                           parse_dates=['Date'])
        except Exception as err:
            raise ValueError('Unable to read Met extracted data file {}. {}'.format(in_file, err))

        if self.verbose > 1:
            print('Hourly dataframe: \n {}'.format(hourly_dataframe))
            print('Hourly dataframe data types: \n {}'.format(hourly_dataframe.dtypes))

        if impute_data:
            # set the imputer options (if we are using them)
            self.imputer = IterativeImputer(random_state=random_state, add_indicator=add_indicator,
                                            initial_strategy=initial_strategy, max_iter=max_iter, verbose=self.verbose,
                                            estimator=estimator)
        self.impute_data = impute_data

        # set the power transform options
        self.transformer = preprocessing.PowerTransformer(method=transformer_method, standardize=transformer_standardize)

        # pull out the daily mean and max values for the site list
        # postprocessing the data set, to get daily data
        daily_dataframe = self.postprocess_organisation(hourly_dataframe)
        # sort the data
        daily_dataframe = daily_dataframe.sort_index()

        if save_to_csv:
            # write this dataset to file
            daily_dataframe.to_csv(self.file_out, index=True, header=True, float_format='%.2f')

        return daily_dataframe

    def load_emep_data(self, filename):
        # load the EMEP model data, or create an empty dataframe (required for logic checks in the workflow)
        if filename is not None:
            filename = Path(filename)
            print('reading emep file')
            try:
                emep_dataframe = pd.read_csv(filename.name)
            except Exception as err:
                raise ValueError('Error loading the emap data from filename: {} . {}'.format(filename, err))
            try:
                return emep_dataframe.rename(columns={'NOx': 'NOXasNO2'})
            except Exception as err:
                raise ValueError('EMEP file does not contain an \'NOx\' column')
        else:
            return pd.DataFrame()

    # functions for station indentifying

    def return_closest_station(self, station_distances):

        stat_info = station_distances.idxmin(axis=0)
        station_data = station_distances.loc[stat_info[0]]

        return station_data.name, station_data[0]

    def get_closest_station_data(self, stations_centre, station_distances, station_id, station_number):

        for ii in range(0, station_number):
            stat_string = 'Station' + str(ii + 1)
            dist_string = 'Distance' + str(ii + 1)

            new_info = self.return_closest_station(station_distances)
            if (new_info[1] == 0):
                station_distances = station_distances.drop(index=new_info[0])
                new_info = self.return_closest_station(station_distances=station_distances)

            stations_centre.loc[station_id, stat_string] = new_info[0]
            stations_centre.loc[station_id, dist_string] = new_info[1]

            station_distances = station_distances.drop(index=new_info[0])

        return stations_centre

    def proc_all_stations(self, stations_centre, stations_input, station_number):

        for ii in range(0, station_number):
            stat_string = 'Station' + str(ii + 1)
            dist_string = 'Distance' + str(ii + 1)
            stations_centre[stat_string] = ''
            stations_centre[dist_string] = np.nan

        for index, row in stations_centre.iterrows():
            print('processing station ', index)
            stat_location = (row['Latitude'], row['Longitude'])
            station_distances = self.calc_station_distances(stations_in=stations_input, stat_location=stat_location)
            stations_centre = self.get_closest_station_data(stations_centre=stations_centre, \
                                                       station_distances=station_distances, \
                                                       station_id=index, \
                                                       station_number=station_number)

        return stations_centre

    def station_listing(self, grouped_data_in):
        '''
        arguments:
            grouped_data_in:
                measurement site data grouped by 24 hour period
                this will be pre-filtered according to what level of
                data you want to keep (so it will be a count of the
                days which meet that criteria)
            min_years (default 1):
                minimum number of years of data that a site must have
            min_years_reference (default 3.5):
                minimum number of years of data for any site that we
                are going to use as a reference site later

        returns:
            required_site_list:
                list of sites with a data count > min_years
            useful_site_list:
                list of sites with a data count > min_years_reference
        '''

        site_list_interior = grouped_data_in.index.levels[0]

        required_site_list = []
        useful_site_list = []

        for site in site_list_interior:
            try:
                date_num = len(grouped_data_in.loc[(site,),])
            except:
                date_num = 0
            if date_num > self.min_years * 365:
                required_site_list.append(site)
                print('\t{} has {} years of data'.format(site, date_num / 365))
            if date_num > self.min_years_reference * 365:
                useful_site_list.append(site)

        return required_site_list, useful_site_list

    def get_station_distances(self, stations_in, site_in, useful_sites_in):

        station_location = stations_in.loc[site_in]['Latitude'], stations_in.loc[site_in]['Longitude']
        station_distances = self.calc_station_distances(stations_in=stations_in.loc[useful_sites_in],
                                                   stat_location=station_location)

        # sort by distance, then drop any station which is the same location as our site of interest
        station_distances = station_distances.sort_values(by='Distance', ascending=True)
        station_distances[station_distances.Distance == 0] = np.nan
        station_distances = station_distances.dropna()

        return station_distances

    def postprocess_data(self, input_dataframe, site):

        working_dataframe = input_dataframe.drop(columns=AurnModule.SITE_ID_NEW)
        tempgroups = working_dataframe.groupby(pd.Grouper(key='Date', freq='1D'))

        data_counts = tempgroups.count()
        data_max = tempgroups.max()
        data_mean = tempgroups.mean()

        cols_old = data_counts.columns

        cols_counts = dict((key, key + '_count') for key in cols_old.values)
        cols_max = dict((key, key + '_max') for key in cols_old.values)
        cols_mean = dict((key, key + '_mean') for key in cols_old.values)

        data_counts = data_counts.rename(columns=cols_counts)
        data_max = data_max.rename(columns=cols_max)
        data_mean = data_mean.rename(columns=cols_mean)

        data_out = data_mean.join([data_max, data_counts])

        # add the site as a new column, and set as part of multiindex with the date
        site_name = "{} [AQ]".format(site)

        data_out['SiteID'] = site_name
        data_out = data_out.reset_index(drop=False).set_index(['Date', 'SiteID'])

        return data_out

    #  testing the reshaping code
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

        if self.verbose > 2: print('df_work input to power transformer: \n {}'.format(df_work))
        # power transformer fitting and transforming
        self.transformer.fit(df_work.dropna())
        if self.verbose > 2: print('Power transformer: Completed data fitting. Beginning power transformation')
        np_out = self.transformer.transform(df_work)
        if self.verbose > 2: print('Power transformer: Completed transformation. Beginning imputation')

        # impute the missing values in this new dataframe
        self.imputer.fit(np_out)
        if self.verbose > 2: print('Imputer: Completed imputation fitting. Beginning imputer tranformation')
        imp_out = self.imputer.transform(np_out)
        if self.verbose > 2: print('Imputer Completed transformation. Beginning inverse transformation')

        # apply the inverse transformation for our datasets (leaving out the indicator flags)
        np_inv = self.transformer.inverse_transform(imp_out[:, :np_out.shape[1]])
        if self.verbose > 2: print('Imputer Completed inverse transformation. Beginning copying and tranforming values')

        # copy the transformed values to a new dataframe
        df_out = df_in.copy(deep=True)
        for pos, col in enumerate(col_save):
            pos_out = list(cols).index(col)
            df_out.iloc[:, pos_out] = np_inv[:, pos]
        if self.verbose > 1: print('Imputation: copied transformed values into new dataframe')

        return df_out

    def postprocess_organisation(self, hourly_dataframe):

        final_dataframe = pd.DataFrame()
        site_list_internal = hourly_dataframe["SiteID"].unique()

        date_index = pd.date_range(start=self.start, end=self.end, freq='1H')
        hourly_dataframe_internal = hourly_dataframe.set_index('Date')

        # do some analysis of the data, getting data counts
        tempgroups = hourly_dataframe.groupby(['SiteID', pd.Grouper(key='Date', freq='1D')])
        daily_hour_counts = tempgroups.count()
        spc_list = daily_hour_counts.columns.values

        # imputation of the values requires more preprocessing and work...
        if self.impute_data:
            # Set station number
            station_number = min(5, len(site_list_internal) - 1)

            req_sites = {}
            use_sites = {}

            for spc in spc_list:
                print('site day counts for {}'.format(spc))
                req_days_counts = daily_hour_counts[spc]
                req_days_counts = req_days_counts[req_days_counts > 0]
                req_sites[spc], use_sites[spc] = self.station_listing(req_days_counts)
                print('VERBOSE: ', self.verbose)
                if self.verbose > 0: print('\t\treq sites {}:'.format(spc), req_sites[spc])
                if self.verbose > 0: print('\t\tuse sites {}:'.format(spc), use_sites[spc])

            if not self._emep_data.empty:
                if self.verbose > 0: print('Loading EMEP data')
                emep_dataframe_internal = self._emep_data.set_index('Date')

            if self.verbose > 1: print('1. Site list internal: ', site_list_internal)
            for site in site_list_internal:
                if self.verbose > 1: print('2. Site: ', site)

                # get list of chemical species that we need to impute for this site (including Date info)
                req_spc = []
                for spc in spc_list:
                    if site in req_sites[spc]:
                        req_spc.append(spc)

                # copy these to a new dataframe
                working_hourly_dataframe = pd.DataFrame([], index=date_index)
                working_hourly_dataframe[req_spc] = \
                    hourly_dataframe_internal[hourly_dataframe_internal['SiteID'] == site][req_spc]

                # get list of neighbouring sites for each of the chemical species of interest
                for spc in spc_list:
                    if self.verbose > 1: print('3. Species: ', spc)
                    station_distances = self.get_station_distances(self.station_data, site, use_sites[spc])
                    if self.verbose > 1: print('4. Station number:', station_number)
                    if self.verbose > 1: print('5. distances:', station_distances)
                    if self.verbose > 1: print('6.', len(station_distances))
                    for ii in range(0, min(station_number, len(station_distances))):
                        if self.verbose > 1: print('7. ii', ii)
                        station_code = station_distances.index[ii]
                        working_hourly_dataframe['{}_{}'.format(spc, station_code)] = \
                            hourly_dataframe_internal[hourly_dataframe_internal['SiteID'] == station_code][spc]

                # get EMEP predictions of chemical species of interest (if needed)
                if self.verbose > 1: print('EMEP data: {}'.format(self._emep_data))
                if not self._emep_data.empty:
                    if self.verbose > 0: print('Using EMEP data')
                    for spc in spc_list:
                        working_hourly_dataframe['{}_{}'.format(spc, 'EMEP')] = \
                            emep_dataframe_internal[emep_dataframe_internal['SiteID'] == site][spc]

                # run the imputation process
                imputed_hourly_dataframe = self.transform_and_impute_data(working_hourly_dataframe)

                # copy imputed data of interest into original dataframe
                for spc in spc_list:
                    working_hourly_dataframe['imputed {}'.format(spc)] = 0
                    if spc in req_spc:
                        working_hourly_dataframe['imputed {}'.format(spc)] = working_hourly_dataframe[spc].isna() * 1
                        working_hourly_dataframe[spc] = imputed_hourly_dataframe[spc]
                    else:
                        working_hourly_dataframe[spc] = np.nan

                # postprocess the new datasets, to get daily mean and max, and populate final dataframe
                daily_grouped_data = working_hourly_dataframe.groupby(pd.Grouper(freq='1D'))
                temp_dataframe = pd.DataFrame()
                for spc in spc_list:
                    temp_dataframe['{}.mean'.format(spc)] = daily_grouped_data.mean()[spc]
                    temp_dataframe['{}.max'.format(spc)] = daily_grouped_data.max()[spc]
                    if spc in req_spc:
                        temp_dataframe['{}.flag'.format(spc)] = daily_grouped_data.mean()['imputed {}'.format(spc)]
                    else:
                        temp_dataframe['{}.flag'.format(spc)] = 0.0
                    temp_dataframe['sensor_name'] = '{} [AQ]'.format(site)

                temp_dataframe = temp_dataframe.rename_axis('time_stamp', axis=0).set_index('sensor_name', append=True)
                final_dataframe = final_dataframe.append(temp_dataframe)

        else:  # simpler post processing of data
            for site in site_list_internal:
                # select our subset of metadata for this station
                station_name = self.station_data.loc[site]['site_name']
                print("processing site {} ({})".format(site, station_name))

                working_hourly_dataframe = pd.DataFrame([], index=date_index)
                working_hourly_dataframe = hourly_dataframe_internal[hourly_dataframe_internal["SiteID"] == site]

                # postprocessing the data set, to get daily data
                daily_grouped_data = working_hourly_dataframe.groupby(pd.Grouper(freq='1D'))
                temp_dataframe = pd.DataFrame()
                for spc in spc_list:
                    temp_dataframe['{}.mean'.format(spc)] = daily_grouped_data.mean()[spc]
                    temp_dataframe['{}.max'.format(spc)] = daily_grouped_data.max()[spc]
                    temp_dataframe['{}.flag'.format(spc)] = 0.0
                    temp_dataframe['sensor_name'] = '{} [AQ]'.format(site)

                temp_dataframe = temp_dataframe.rename_axis('time_stamp', axis=0).set_index('sensor_name', append=True)
                final_dataframe = final_dataframe.append(temp_dataframe)

        return final_dataframe


    #  testing the reshaping code

    def test_preprocess_code(self, df_in, spc_zero_process=['O3', 'NO2', 'NOXasNO2'], min_value=0.01):

        # define the species which we will convert <=0 values to a minimum value
        # spc_zero_process = ['O3','NO2','NOXasNO2']
        # min_value = 0.01

        # copy the input array, and note the columns
        df_work = df_in.copy(deep=True)
        cols = df_in.columns

        # find missing datasets to remove, and clean up negative and zero values
        # also we note the columns that will be saved, for transferring data back!
        col_remove = []
        col_save = []
        for col in cols:
            if all(df_work[col].isna()):
                col_remove.append(col)
            else:
                col_save.append(col)
                if (col in spc_zero_process):
                    df_work[col][df_work[col] <= 0.0] = min_value
                else:
                    df_work[col][df_work[col] <= 0.0] = np.nan
        df_work = df_work.drop(columns=col_remove)

        # power transformer fitting and transforming
        self.transformer.fit(df_work.dropna())
        np_out = self.transformer.transform(df_work)

        # copy the transformed values to a new dataframe
        df_out = df_in.copy(deep=True)
        for pos, col in enumerate(col_save):
            pos_out = list(cols).index(col)
            df_out.iloc[:, pos_out] = np_out[:, pos]

        return self.transformer, df_out
