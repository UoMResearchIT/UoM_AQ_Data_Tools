try:
    from datetime import datetime
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
    """
        Class used for post-processing data that has been extracted from AURN server.
    """
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
        """ Initialise instance of the AurnPostProcessor class.
            Initialises the private class variables

            Args:
                metadata_filename: filename of the metadata used in Aurn data extraction
                metadata_url: alternative source of AURN metadata, if metadata_filename is None
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of AurnPostProcessor

        """
        super(AurnPostProcessor, self).__init__(out_dir, verbose)
        AurnModule.__init__(self, metadata_filename=metadata_filename, metadata_url=metadata_url)
        DateRangeProcessor.__init__(self)

        self._emep_data = None
        self.min_years_reference = AurnPostProcessor.DEFAULT_MIN_YEARS_REFERENCE
        self.min_years = AurnPostProcessor.DEFAULT_MIN_YEARS
        self.impute_data = False
        self._imputer = None
        self._transformer = None

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
            station_data = station_data.set_index('site_id')
        except Exception as err:
            raise ValueError('Unable to get correct site data from Metadata input file. Check metadata file content.')

        self._station_data = station_data


    def impute_method_setup(self, random_state=DEFAULT_IMPUTER_RANDOM_STATE, add_indicator=DEFAULT_IMPUTER_ADD_INDICATOR,
                 initial_strategy=DEFAULT_IMPUTER_INITIAL_STRATEGY,
                 max_iter=DEFAULT_IMPUTER_MAX_ITER, estimator=DEFAULT_IMPUTER_ESTIMATOR,
                 transformer_method=DEFAULT_TRANSFORMER_METHOD, transformer_standardize=DEFAULT_TRANSFORMER_STANDARDIZE):
        """ Initialises the IterativeImputer and PowerTransformer methods required if missing data is to be imputed.
            Parameters are passed to the sklearn routines. For further documentation on how these functions work, 
            and what the parameters denote, please refer to the sklearn documentation.

            IterativeImputer: https://scikit-learn.org/stable/modules/generated/sklearn.impute.IterativeImputer.html
            PowerTransformer: https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.PowerTransformer.html
            
            Args:
                random_state:           (int) (IterativeImputer) seed for pseudo random number generator
                add_indicator:          (boolean) (IterativeImputer) if True adds a `MissingIndicator` transform to the stack
                initial_strategy:       (str) (IterativeImputer) define strategy to use for initialising missing values
                max_iter:               (int) (IterativeImputer) maximum number of imputation rounds to perform
                estimator:              (str) (IterativeImputer) estimator method to be used
                transformer_standardize:(boolean) (PowerTransformer) apply zero-mean, unit-variance normalization
                transformer_method:     (str) (PowerTransformer) power transform method to use

             Returns: None
        """

        # set the imputer options (if we are using them)
        self.imputer = IterativeImputer(random_state=random_state, add_indicator=add_indicator,
                                        initial_strategy=initial_strategy, max_iter=max_iter, verbose=self.verbose,
                                        estimator=estimator)

        # set the power transform options
        self.transformer = preprocessing.PowerTransformer(method=transformer_method, standardize=transformer_standardize)




    def process(self, in_file, date_range=None,
                site_list=AurnModule.DEFAULT_SITE_LIST,
                emep_filename=DEFAULT_EMEP_FILENAME,
                min_years_reference=DEFAULT_MIN_YEARS_REFERENCE,
                min_years=DEFAULT_MIN_YEARS,
                impute_data=PostProcessor.DEFAULT_IMPUTE_DATA,
                save_to_csv=PostProcessor.DEFAULT_SAVE_TO_CSV,
                outfile_suffix=''):

        """ Post process the data extracted from the AURN dataset, based on the parameters given.
            
            Args:
                in_file:                (str) The file spec of the input file (required)
                date_range:             (list of 2 datetime) The date range of interest
                site_list:              (list of string/number) Site IDs of interest
                emep_filename:          (str) The file spec of the EMEP file to be used to help calculate #Todo Doug
                min_years_reference:    (float) The minimum number of years of data for any site that  we are going to
                                            use as a reference site later. (this cannot be less than min_years)
                min_years:              (float) The minimum number of years of data that a site must have
                impute_data:            (boolean) Whether to attempt to impute missing data
                save_to_csv:            (boolean) Whether to save the output dateframes to CSV file(s)
                outfile_suffix:         (str) The suffix to appended to the end of output file names.

            Returns:
                daily_dataframe: daily dataset, for all measurements, as pandas.Dataframe
                    Required MultiIndex:
                        'time_stamp'  (datetime object): date (only) (e.g. 2017-06-01)
                        'sensor_name'          (string): ID string for site (e.g. 'LIN3 [AQ]')
                    Required columns:
                        'O3.max'       (float): daily maximum value
                        'O3.mean'      (float): daily mean value
                        'O3.flag'      (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
                        'PM10.max'       (float): daily maximum value
                        'PM10.mean'      (float): daily mean value
                        'PM10.flag'      (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
                        'PM2.5.max'       (float): daily maximum value
                        'PM2.5.mean'      (float): daily mean value
                        'PM2.5.flag'      (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
                        'NO2.max'       (float): daily maximum value
                        'NO2.mean'      (float): daily mean value
                        'NO2.flag'      (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
                        'NOXasNO2.max'       (float): daily maximum value
                        'NOXasNO2.mean'      (float): daily mean value
                        'NOXasNO2.flag'      (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
                        'SO2.max'       (float): daily maximum value
                        'SO2.mean'      (float): daily mean value
                        'SO2.flag'      (float): flag to indicate fraction of imputed data
                                                        (1 = fully imputed, 0 = no imputed values were used)
        """

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

        self.impute_data = impute_data

        # load and prepare the hourly dataset
        hourly_dataframe = self.load_aurn_data(in_file)

        print('filter for minimum data lengths, and reduce dataset to only stations of interest')
        hourly_dataframe_filtered, reference_sites, required_sites, site_list_internal = \
            self.list_required_and_reference_sites(hourly_dataframe)
        # get the list of required sites from what is available, and what was requested
        site_list_internal = set(site_list_internal).intersection(self.site_list)

        if len(hourly_dataframe_filtered.index) == 0:
            print('Exiting post-processing: Metadata is empty after initial filtering processes')
            return

        if self.impute_data:
            print('imputation of data, returning hourly data')
            hourly_dataframe = self.organise_data_imputation(
                hourly_dataframe_filtered, reference_sites, required_sites, site_list_internal)
        else:
            print('sorting data (no imputation), returning hourly data')
            hourly_dataframe = self.organise_data(hourly_dataframe_filtered, site_list_internal)

        # calculate the daily max and mean for each station
        daily_dataframe = self.combine_and_organise_mean_max(hourly_dataframe)

        if save_to_csv:
            # write this dataset to file
            daily_dataframe.to_csv(self.file_out, index=True, header=True, float_format='%.2f')

        return daily_dataframe

    def combine_and_organise_mean_max(self, hourly_dataframe):
        """
        Combine and organise the daily mean, maximum, and count information.
        
        Args:
            hourly_dataframe: hourly dataset, for all measurements, as pandas.Dataframe
                Required Index:
                    Date   (datetime object):
                    SiteID          (string):
                Optional Columns:
                    O3       (float):
                    PM10     (float):
                    PM2.5    (float):
                    NO2      (float):
                    NOXasNO2 (float):
                    SO2      (float):
                    imputed O3       (int): flag indicating imputed data (0=original,1=imputed)
                    imputed PM10     (int):
                    imputed PM2.5    (int):
                    imputed NO2      (int):
                    imputed NOXasNO2 (int):
                    imputed SO2      (int):
            
        Returns:
            final_dataframe: daily dataset, for all measurements, as pandas.Dataframe
                Required MultiIndex:
                    'time_stamp'  (datetime object): date (only) (e.g. 2017-06-01)
                    'sensor_name'          (string): ID string for site (e.g. 'LIN3 [AQ]')
                Required columns:
                    'O3.max'       (float): daily maximum value
                    'O3.mean'      (float): daily mean value
                    'O3.flag'      (float): flag to indicate fraction of imputed data
                                                    (1 = fully imputed, 0 = no imputed values were used)
                    'PM10.max'       (float): daily maximum value
                    'PM10.mean'      (float): daily mean value
                    'PM10.flag'      (float): flag to indicate fraction of imputed data
                                                    (1 = fully imputed, 0 = no imputed values were used)
                    'PM2.5.max'       (float): daily maximum value
                    'PM2.5.mean'      (float): daily mean value
                    'PM2.5.flag'      (float): flag to indicate fraction of imputed data
                                                    (1 = fully imputed, 0 = no imputed values were used)
                    'NO2.max'       (float): daily maximum value
                    'NO2.mean'      (float): daily mean value
                    'NO2.flag'      (float): flag to indicate fraction of imputed data
                                                    (1 = fully imputed, 0 = no imputed values were used)
                    'NOXasNO2.max'       (float): daily maximum value
                    'NOXasNO2.mean'      (float): daily mean value
                    'NOXasNO2.flag'      (float): flag to indicate fraction of imputed data
                                                    (1 = fully imputed, 0 = no imputed values were used)
                    'SO2.max'       (float): daily maximum value
                    'SO2.mean'      (float): daily mean value
                    'SO2.flag'      (float): flag to indicate fraction of imputed data
                                                    (1 = fully imputed, 0 = no imputed values were used)
        """


        #### group by date and site
        daily_grouped_data = hourly_dataframe.groupby([pd.Grouper(level="Date", freq='1D'), 'SiteID'])
        spc_list = ['O3', 'PM10', 'PM2.5', 'NO2', 'NOXasNO2', 'SO2'] # TODO Doug - make this check database columns! 
        
        
        #### loop by spc through grouped data, and calculate the mean, max, and flag values
        for spc in spc_list:
            temp_dataframe = pd.DataFrame()
            temp_dataframe['{}_mean'.format(spc)] = daily_grouped_data.mean()[spc]
            temp_dataframe['{}_max'.format(spc)] = daily_grouped_data.max()[spc]
            temp_dataframe['{}_flag'.format(spc)] = daily_grouped_data.mean()['{}_flag'.format(spc)]
            try:
                final_dataframe = final_dataframe.merge(temp_dataframe, how='outer', left_index=True, right_index=True)
            except:
                final_dataframe = temp_dataframe.copy()

        #### rename the sites, to include AQ flag
        final_dataframe.index = final_dataframe.index.set_levels(
                    ['{} [AQ]'.format(x) for x in final_dataframe.index.levels[1]], level=1)

        #### rename the index columns
        final_dataframe.index.rename(['time_stamp', 'sensor_name'], inplace=True)

        #### return output dataframe
        return(final_dataframe)

    def load_aurn_data(self, file_in):
        """
        Loading the AURN dataset.
        
        Args:
            file_in (Path object or string): path for the file to be read in
            
        Returns:
            hourly_dataframe:  hourly dataset, for all measurements, as pandas.Dataframe
                    Index: none
                    Required Columns:
                        Date   (datetime object):
                        SiteID          (string):
                    Optional Columns:
                        O3       (float):
                        PM10     (float):
                        PM2.5    (float):
                        NO2      (float):
                        NOXasNO2 (float):
                        SO2      (float):
        """
        # Read in hourly dataframe file
        try:
            hourly_dataframe = pd.read_csv(file_in,
                                           sep=',',
                                           usecols=[AurnModule.INDEX_EXTRACTED].append(AurnModule.NEW_FILE_COLS),
                                           index_col=AurnModule.INDEX_EXTRACTED,
                                           parse_dates=['Date'])
        except Exception as err:
            raise ValueError('Unable to read Met extracted data file {}. {}'.format(file_in, err))

        if self.verbose > 1:
            print('Hourly dataframe: \n {}'.format(hourly_dataframe))
            print('Hourly dataframe data types: \n {}'.format(hourly_dataframe.dtypes))
            
        return(hourly_dataframe)

    def load_emep_data(self, filename):
        """
        Loads the EMEP model data, or create an empty dataframe (required for logic checks in the workflow)
        
        Args:
            filename (str): location of the EMEP file. This should be empty if there is no EMEP data
        
        Returns:
            emep_dataframe: pandas Dataframe, containing the EMEP model data. If no EMEP data
                            is to be used then this will be an empty Dataframe.
                    Index: none
                    Required Columns:
                        Date   (datetime object):
                        SiteID          (string):
                        O3       (float):
                        PM10     (float):
                        PM2.5    (float):
                        NO2      (float):
                        NOXasNO2 (float):
                        SO2      (float):
        """
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

    def list_required_and_reference_sites(self, data_in):
        """
        This function creates the lists of required sites, and reference sites, for the 
        final dataset.
        
        Args:
            data_in: hourly dataset, for all measurements, as pandas.Dataframe
                Index: none
                Required Columns:
                    Date   (datetime object):
                    SiteID          (string):
                Optional Columns:
                    O3       (float):
                    PM10     (float):
                    PM2.5    (float):
                    NO2      (float):
                    NOXasNO2 (float):
                    SO2      (float):
            
        Returns:
            met_data_filtered: pandas dataframe, as above, containing hourly dataset for only 
                               the required station datasets
            reference_sites: (dict, keys are species):
                            items: (list of strings) the siteID's for our reference sites for each `spc` 
            required_sites: (dict, keys are species):
                            items: (list of strings) required sites for `spc`
            combined_req_site_list: (list, strings) a single list of required sites
        
        """
        print('    get the lists of required and reference stations for each measurement variable')
        tempgroups = data_in.groupby(['SiteID', pd.Grouper(key='Date', freq='1D')])
        daily_hour_counts = tempgroups.count()
        spc_list = daily_hour_counts.columns.values
        
        required_sites = {}
        reference_sites = {}
        combined_req_site_list = []
        
        for spc in spc_list:
            print('site day counts for {}'.format(spc))
            req_days_counts = daily_hour_counts[spc]
            req_days_counts = req_days_counts[req_days_counts > 0]
            required_sites[spc], reference_sites[spc] = self.station_listing(req_days_counts)
            combined_req_site_list = combined_req_site_list + required_sites[spc]

            print('VERBOSE: ', self.verbose)
            if self.verbose > 0: print('\t\treq sites {}:'.format(spc), required_sites[spc])
            if self.verbose > 0: print('\t\tuse sites {}:'.format(spc), reference_sites[spc])

        # get a list of all sites which are required for at least one measurement set
        combined_req_site_list = list(dict.fromkeys(combined_req_site_list))
        print('there are {} required sites, and {} reference sites'.format(len(required_sites), len(reference_sites)))
        data_filtered = data_in[data_in['SiteID'].isin(combined_req_site_list)]

        return data_filtered, reference_sites, required_sites, combined_req_site_list

    def organise_data_imputation(self, hourly_dataframe_filtered, reference_sites, required_sites, site_list_internal): 
        """
        Function for organising the imputation of the datasets. This runs the 
        'transform_and_impute_data' function for each of the variables of interest.
        
        Args:
            hourly_dataframe_filtered: hourly dataset, for all measurements, as pandas.Dataframe
                Index: none
                Required Columns:
                    Date   (datetime object):
                    SiteID          (string):
                Optional Columns:
                    O3       (float):
                    PM10     (float):
                    PM2.5    (float):
                    NO2      (float):
                    NOXasNO2 (float):
                    SO2      (float):
            reference_sites (list, string or int): sites to use for reference when imputing datasets
            required_sites: (dict, keys are species):
                            items: (list of strings) required sites for `spc`
            site_list_internal (list, string or int): combined list of sites to retain
        
        Returns:
            output_dataframe: hourly dataset, for all measurements, as pandas.Dataframe
                Required Index:
                    Date   (datetime object):
                    SiteID          (string):
                Optional Columns:
                    O3       (float):
                    PM10     (float):
                    PM2.5    (float):
                    NO2      (float):
                    NOXasNO2 (float):
                    SO2      (float):
                    O3_flag       (int): flag indicating imputed data (0=original,1=imputed)
                    PM10_flag     (int):
                    PM2.5_flag    (int):
                    NO2_flag      (int):
                    NOXasNO2_flag (int):
                    SO2_flag      (int):
        """

        output_dataframe = pd.DataFrame()
        date_index = pd.date_range(start=self.start, end=self.end, freq='1H', name='Date')

        # Set the number of reference stations to request
        ref_station_numbers = [len(reference_sites[x]) for x in reference_sites.keys()] 
        station_number = min([5] + [len(ref_station_numbers) - 1])
        
        hourly_dataframe_internal = hourly_dataframe_filtered.set_index('Date')
        spc_list = ['O3','PM10','PM2.5','NO2','NOXasNO2','SO2'] # TODO Doug - make this check database columns! 


        if not self._emep_data.empty:
            if self.verbose > 0: print('Loading EMEP data')
            emep_dataframe_internal = self._emep_data.set_index('Date')

        if self.verbose > 1: print('1. Site list internal: ', site_list_internal)
        for site in site_list_internal:
            if self.verbose > 1: print('2. Site: ', site)

            # get list of chemical species that we need to impute for this site (including Date info)
            req_spc = []
            for spc in spc_list:
                if site in required_sites[spc]:
                    req_spc.append(spc)

            # copy these to a new dataframe
            working_hourly_dataframe = pd.DataFrame([], index=date_index)
            working_hourly_dataframe[req_spc] = \
                hourly_dataframe_internal[hourly_dataframe_internal['SiteID'] == site][req_spc]
            copy_hourly_dataframe = working_hourly_dataframe.copy()
            copy_hourly_dataframe['SiteID'] = site

            # get list of neighbouring sites for each of the chemical species of interest
            for spc in spc_list:
                if self.verbose > 1: print('3. Species: ', spc)
                station_distances = self.get_station_distances(site, reference_sites[spc])
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

            # copy imputed data of interest into copy of original dataframe (without EMEP and neighbouring sites)
            for spc in spc_list:
                copy_hourly_dataframe['{}_flag'.format(spc)] = 0
                if spc in req_spc:
                    copy_hourly_dataframe['{}_flag'.format(spc)] = copy_hourly_dataframe[spc].isna() * 1
                    copy_hourly_dataframe[spc] = imputed_hourly_dataframe[spc]
                else:
                    copy_hourly_dataframe[spc] = np.nan

            output_dataframe = output_dataframe.append(copy_hourly_dataframe)

        output_dataframe = output_dataframe.reset_index().set_index(['Date','SiteID'])
        return(output_dataframe)

    def organise_data(self, hourly_dataframe_filtered, site_list_internal): 
        """
        Function for organising the required datasets. This mirrors the imputation function.
        
        Args:
            hourly_dataframe_filtered: hourly dataset, for all measurements, as pandas.Dataframe
                Index: none
                Required Columns:
                    Date   (datetime object):
                    SiteID          (string):
                Optional Columns:
                    O3       (float):
                    PM10     (float):
                    PM2.5    (float):
                    NO2      (float):
                    NOXasNO2 (float):
                    SO2      (float):
            site_list_internal (list, string or int): combined list of sites to retain
        
        Returns:
            hourly_dataframe: hourly dataset, for all measurements, as pandas.Dataframe
                Required Index:
                    Date   (datetime object):
                    SiteID          (string):
                Optional Columns:
                    O3       (float):
                    PM10     (float):
                    PM2.5    (float):
                    NO2      (float):
                    NOXasNO2 (float):
                    SO2      (float):
                    O3_flag       (int): flag indicating imputed data (0=original,1=imputed)
                    PM10_flag     (int):
                    PM2.5_flag    (int):
                    NO2_flag      (int):
                    NOXasNO2_flag (int):
                    SO2_flag      (int):
        """

        date_index = pd.date_range(start=self.start, end=self.end, freq='1H', name='Date')
        output_dataframe = pd.DataFrame()

        hourly_dataframe_internal = hourly_dataframe_filtered.set_index('Date')
        spc_list = ['O3','PM10','PM2.5','NO2','NOXasNO2','SO2'] # TODO Doug - make this check database columns! 

        if self.verbose > 1: print('1. Site list internal: ', site_list_internal)
        for site in site_list_internal:
            if self.verbose > 1: print('2. Site: ', site)

            # create new dataframe, with the dates that we are interested in
            working_hourly_dataframe = pd.DataFrame([], index=date_index)
            working_hourly_dataframe['SiteID'] = site

            # copy these to a new dataframe
            working_hourly_dataframe[spc_list] = \
                hourly_dataframe_internal[hourly_dataframe_internal['SiteID'] == site][spc_list]

            # copy imputed data of interest into copy of original dataframe (without EMEP and neighbouring sites)
            for spc in spc_list:
                working_hourly_dataframe['{}_flag'.format(spc)] = 0

            # append data to the output dataframe
            output_dataframe = output_dataframe.append(working_hourly_dataframe)

        output_dataframe = output_dataframe.reset_index().set_index(['Date','SiteID'])
        return(output_dataframe)

    def transform_and_impute_data(self, df_in):
        """
        Function for organising the transformation of the dataset, then imputing missing
        data, before detransforming the data and returning it.
        
        Args:
            df_in: pandas dataframe containing the datasets to impute
                Required Index:
                    date (datetime64 objects): date / time for each reading
                Optional Columns: Measurement data at the site for which we are imputing
                                  the data. Only those pollutants which have been measured
                                  at this site will be included.
                    O3       (float): 
                    PM10     (float):
                    PM2.5    (float):
                    NO2      (float):
                    NOXasNO2 (float):
                    SO2      (float):
                Reference Columns: Reference data at the X nearest sites to the 
                                   measurement being processed. All datasets will be
                                   included, even for those pollutants which were not
                                   included in the optional columns above. So, if
                                   5 reference stations are used, this will give 30 (5*6)
                                   columns of reference data. If EMEP data is being used
                                   then these are added for EMEP data too, but only at 
                                   the station of interest (so only another 6 columns are
                                   added). 
                    O3_[siteID]       (float): 
                    PM10_[siteID]     (float):
                    PM2.5_[siteID]    (float):
                    NO2_[siteID]      (float):
                    NOXasNO2_[siteID] (float):
                    SO2_[siteID]      (float):
             
        Returns:
            df_out: pandas dataframe, containing the same datasets as above, but including
                    the imputed data too. All imputed data is included (including that for
                    the reference sites) - it is the task of the calling function to only
                    retain the imputed data for the station of interest, and to discard 
                    the rest of the imputed data.
        """
        
        
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
