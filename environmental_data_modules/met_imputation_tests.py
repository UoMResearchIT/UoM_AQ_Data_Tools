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
    
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    
    from scipy import stats
    
    sns.set_theme()

except:
    pass  # print('Warning: Unable to load library: {}'.format(err))

from environmental_data_modules import MetPostProcessor


class MetImputationTest(MetPostProcessor):
    """
        Class used testing the imputation of the data that is extracted from the MEDMI server.
    """

    # testing default values
    DEFAULT_DATA_LOST = 0.5
    DEFAULT_DATA_LOSS_POSITION = 'end'
    DEFAULT_CHECK_SITES = False
    DEFAULT_SITE_LIST = None
    
    BASE_IMPUTED_STATS_PDF_FILE = '{}/{}_{}_imputed_comparison.pdf'
    BASE_IMPUTED_STATS_CSV_FILE = '{}/met_{}_correlation_stats.csv'
    DEFAULT_FLOAT_FORMAT = '%.4f'


    def __init__(self, out_dir=MetPostProcessor.DEFAULT_OUT_DIR, station_data_filename=MetPostProcessor.DEFAULT_STATION_DATA_FILENAME,
                 verbose=MetPostProcessor.DEFAULT_VERBOSE):
        """ Initialise instance of the MetImputationTest class.
            Initialises the private class variables

            Args:
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of MetPostProcessor

        """
        super(MetImputationTest, self).__init__(out_dir, station_data_filename, verbose)
        self._data_lost = MetImputationTest.DEFAULT_DATA_LOST
        self._data_loss_position = MetImputationTest.DEFAULT_DATA_LOSS_POSITION
        self._site_list = MetImputationTest.DEFAULT_SITE_LIST
        self.check_sites = MetImputationTest.DEFAULT_CHECK_SITES
        self.species_list = MetPostProcessor.SPECIES_PROCESS_LIST
        self.rng = np.random.RandomState(0)

    @property
    def data_lost(self):
        return self._data_lost
    
    @data_lost.setter
    def data_lost(self, data_lost):
        if type(data_lost) == float and 0 <= data_lost < 1:
            self._data_lost = data_lost
        else:
            raise Exception('data_lost is {}, but should be float from 0 upto (not including) 1'.format(data_lost))
            
    @property
    def data_loss_position(self):
        return self._data_loss_position
    
    @data_loss_position.setter
    def data_loss_position(self, data_loss_position):
        if type(data_loss_position) == str and data_loss_position in ['start','middle','end','random']:
            self._data_loss_position = data_loss_position
        else:
            raise Exception("data_loss_position is {}, but should be string matching one of these: ['start','middle','end','random']".format(data_loss_position))

    @property
    def site_list(self):
        return self._site_list

    @site_list.setter
    def site_list(self, site_list):
        """ Set the site_list property.
            If site_list is None, will use the site lists from the AURN metadata property
            Otherwise, will check that each site ID in site_list is in the list of all site IDs in the AURN Metadata

            Dependencies:
                self.metadata

            Args:
                site_list: list of site identifiers (strings or numbers)

            Returns:
                None

        """
        all_sites = self.station_data.index.values
        # get list of sites to process extract_site_data
        if site_list is None:
            self._site_list = all_sites
        else:
            try:
                site_list = set(list(site_list))
            except Exception:
                raise TypeError('Site list must be a list. Input: {}'.format(site_list))
            error_sites = set(site_list) - set(all_sites)
            assert len(error_sites) == 0, ValueError(
                "Each site must be contained in available sites (from Met metadata: {}. Error sites: {}".format(
                    all_sites, str(error_sites)))
            self._site_list = site_list


    def imputation_test(self, in_file, outfile_suffix='', date_range=None,
                exclude_site_list=MetPostProcessor.DEFAULT_EXCLUDE_STATION_LIST,
                min_temperature=MetPostProcessor.DEFAULT_MIN_TEMPERATURE, reference_num_stations=MetPostProcessor.DEFAULT_REFERENCE_NUMBER_STATIONS,
                min_years=MetPostProcessor.DEFAULT_MIN_YEARS, min_years_reference=MetPostProcessor.DEFAULT_MIN_YEARS_REFERENCE,
                data_lost=DEFAULT_DATA_LOST, 
                data_loss_position=DEFAULT_DATA_LOSS_POSITION,
                print_stats=MetPostProcessor.DEFAULT_PRINT_STATS,
                save_to_csv=MetPostProcessor.DEFAULT_SAVE_TO_CSV,
                check_sites=DEFAULT_CHECK_SITES,
                site_list=DEFAULT_SITE_LIST):
        """ Testing the imputation methods used for filling in missing data. Replicates the
            methods used in 'process' function, but for the defined sites will remove a
            given amount of data from the dataframe before imputation, then compare the
            imputed data with the original data, to determine accuracy of the method.
            
            The stations to be tested need to meet the requirements for 'reference' sites.
            
            Args:
                in_file:                (str) The file spec of the input file (required)
                date_range:             (list of 2 datetime) The date range of interest
                exclude_site_list:      (list of string/number) Site IDs to be ignored
                min_years_reference:    (float) The minimum number of years of data for any site that  we are going to
                                            use as a reference site later. (this cannot be less than min_years)
                min_years:              (float) The minimum number of years of data that a site must have
                min_temperature:        (float) The minimum temperature to be used (lower are ignored)
                reference_num_stations: (int) The number of stations to be used for imputation
                data_lost:              (float) The fraction of each dataset to remove
                data_loss_position:     (str) where to lose the data (start,middle,end,random)
                print_stats:            (boolean) Whether to printout the calculation statistics
                save_to_csv:            (boolean) Whether to save the output dateframes to CSV file(s)
                outfile_suffix:         (str) The suffix to appended to the end of output file names.
                check_sites:            (boolean) If True then routine will list appropriate stations
                                            to use for the imputation tests, then exit

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
        self.data_lost = data_lost
        self.data_loss_position = data_loss_position
        self.check_sites = check_sites
        if site_list:
            self.site_list = [int(x) for x in site_list]
        self.pdf_file_string = MetImputationTest.BASE_IMPUTED_STATS_PDF_FILE
        self.csv_file_string = MetImputationTest.BASE_IMPUTED_STATS_CSV_FILE
        self.float_format = MetImputationTest.DEFAULT_FLOAT_FORMAT



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
        met_extracted_data, ref_sites, req_sites, site_list_internal = self.site_list_and_preparation(met_extracted_data)
 
        if len(met_extracted_data.index) == 0:
            print('Exiting post-processing: Metadata is empty after initial filtering processes')
            return
 
        # extract species specific required site lists (for imputation we just use the site_list_internal list)
        req_sites_temp = site_list_internal
        req_sites_pres = site_list_internal
        req_sites_dewpoint = site_list_internal
        # find a unified list of useful sites for all our measurements
        reference_sites = [x for x in ref_sites['dewpoint'] if x in ref_sites['pressure']]



        print('data preparation, to create the test dataset for imputation')
        met_impute_test_data = self.data_preparation(met_extracted_data, 
                                               ref_sites, site_list_internal)

        print('imputation of data, returning hourly data')
        met_data_imputed_temp, met_data_imputed_pres, met_data_imputed_dewpoint = self.organise_data_imputation(
            met_impute_test_data, reference_sites, req_sites_temp, req_sites_pres, req_sites_dewpoint)
        print('sorting data (no imputation), returning hourly data')
        met_data_ref_temp, met_data_ref_pres, met_data_ref_dewpoint = self.organise_data(
            met_extracted_data, req_sites_temp, req_sites_pres, req_sites_dewpoint)

        print('calculation of relative humidity from temperature and dew point temperature for imputed data')
        met_data_imputed_rh = self.rh_calculations(met_impute_test_data, met_data_imputed_temp, met_data_imputed_dewpoint)
        print('calculation of relative humidity from temperature and dew point temperature for reference data')
        met_data_ref_rh = self.rh_calculations(met_extracted_data, met_data_ref_temp, met_data_ref_dewpoint)

        # calculate the daily max and mean for each station
        met_data_imputed_daily = self.combine_and_organise_mean_max(
            met_impute_test_data, met_data_imputed_temp, met_data_imputed_pres, met_data_imputed_rh)
        met_data_ref_daily = self.combine_and_organise_mean_max(
            met_extracted_data, met_data_ref_temp, met_data_ref_pres, met_data_ref_rh)

        # combine hourly datasets to match expected dataframes for plotting data
        hourly_imputed_dataframe = self.combine_plotting_dataset(
            met_impute_test_data, met_data_imputed_temp, met_data_imputed_pres, met_data_imputed_rh)
        hourly_reference_dataframe = self.combine_plotting_dataset(
            met_extracted_data, met_data_ref_temp, met_data_ref_pres, met_data_ref_rh)

        # calculate the stats for the hourly and daily data, printing out graphs of these
        self.imputation_hourly_analysis(hourly_imputed_dataframe,hourly_reference_dataframe,site_list_internal)
        self.imputation_daily_analysis(met_data_imputed_daily,met_data_ref_daily,site_list_internal)


        #if save_to_csv:
        #    # write data to file
        #    if self.verbose > 1: print('Writing to file: {}'.format(self.file_out.format(outfile_suffix)))
        #    met_data_daily.to_csv(self.file_out, index=True, header=True, float_format='%.2f')

        return site_list_internal, met_impute_test_data, met_data_imputed_daily

    def site_list_and_preparation(self,hourly_dataframe):
        """
        Wrapper for the list_required_and_reference_sites routine. This will list sites 
        based on the given minimum and reference year requirements, then determine which
        sites for all species of interest fit the reference year requirements. If none do,
        or if 'check_sites' flag is True, then potential sites of use will be listed, and 
        the program exited.
        
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
                               the reference station datasets
                               
            reference_sites_out: (dict, keys are species):
                            items: (list of strings) the siteID's for our reference sites for each `spc`
                                                     minus the sites in the `site_working` list
            required_sites: (dict, keys are species):
                            items: (list of strings) required sites for `spc` (not used later?)
            site_working_list: (list, strings) a single list of required sites for the imputation tests
            
            reference_sites: (list, string or int) the siteID's for our reference sites 
                             (composite list made up from the measurement specific lists below)
            req_sites_temp: (list, string or int) required sites for temperature data
            req_sites_pres: (list, string or int) required sites for pressure data
            req_sites_dewpoint: (list, string or int) required sites for dewpoint temperature data
        """

        df_part_filtered, reference_sites, required_sites, site_list_internal = self.list_required_and_reference_sites(hourly_dataframe)


        # create a list of the sites that we can use for the imputation tests for all requested species
        site_all = site_list_internal.copy()
        combined_reference_site_list = []
        for spc in self.species_list:
            site_all = set(site_all).intersection(reference_sites[spc])
            combined_reference_site_list = combined_reference_site_list + required_sites[spc]

        # trim down the database to cover only stations which are references for at least one species
        combined_reference_site_list = list(dict.fromkeys(combined_reference_site_list))
        hourly_dataframe_filtered = df_part_filtered[df_part_filtered[self._site_string].isin(combined_reference_site_list)]

        # get the list of required sites from what is available, and what was requested
        site_working_list = set(site_all).intersection(self.site_list)

        # checks on what sites are available, and if we use them or not
        if self.check_sites or len(site_working_list) == 0:
            for spc in self.species_list:
                print('for species {} there are {} sites suitable for testing imputation'.
                         format(spc,len(reference_sites[spc])))
                print(reference_sites[spc])
            print('there are {} sites suitable for all requested species'.format(len(site_all)))
            print(site_all)
            if not self.check_sites:
                print('Requested sites were: {}'.format(self.site_list))
                print('We are exiting because none were suitable sites for imputation tests for all requested species, see above messages.')
            return pd.DataFrame(), [], [], []

        # derive new lists of reference stations, excluding the sites we will use for imputation tests
        reference_sites_out = {}
        for spc in self.species_list:
            reference_sites_out[spc] = [site for site in reference_sites[spc] if site not in site_working_list]

        # success! return the filtered dataframe, and our lists of sites
        return hourly_dataframe_filtered, reference_sites_out, required_sites, site_working_list

    def data_preparation(self, hourly_dataframe, reference_sites, site_list_internal):
        """
        Prepare test data for imputation, by removing the specified amount of data from the test sites.
        
        Args:
            hourly_dataframe: hourly dataset, for all measurements, as pandas.Dataframe
                Index: none
                Required columns:
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
            reference_sites: (dict, keys are species):
                            items: (list of strings) the siteID's for our reference sites for each `spc`
                                                     minus the sites in the `site_working` list
            site_list_internal: (list, strings) a single list of required sites
        
        Returns:
            hourly_dataframe_out: hourly dataset, for all measurements, as pandas.Dataframe, with the required data removed
                Index: none
                    'date'        (datetime object) date/time of measurement
                    'siteID'      (string) ID string for site
                    'temperature' (float): temperature
                    'rel_hum'     (float): relative humidity
                    'pressure'    (float): pressure
                    'dewpoint'    (float): dewpoint temperature
        
        """
        
        # get list of reference sites, which should exclude the sites for testing
        reference_site_list = []
        for spc in self.species_list:
            reference_site_list = reference_site_list + reference_sites[spc]
        reference_site_list = list(dict.fromkeys(reference_site_list))
        
        # create dataframe with reference sites only
        hourly_dataframe_out = hourly_dataframe[hourly_dataframe[self._site_string].isin(reference_site_list)]
        
        for site in site_list_internal:
            print('  filtering site {}'.format(site))
            working_dataframe = hourly_dataframe[hourly_dataframe[self._site_string]==site].copy()
            data_length = len(working_dataframe)
            print('index length is {}'.format(data_length))
            if self.data_loss_position == 'end':
                start_point = 0
                end_point = int(np.floor(data_length * self.data_lost))
                working_dataframe = working_dataframe.iloc[start_point:end_point]
            elif self.data_loss_position == 'middle':
                half_data_retain = (1-self.data_lost)/2
                start_point = int(np.floor(data_length * half_data_retain))
                end_point = data_length - start_point
                working_dataframe_start = working_dataframe.iloc[0:start_point]
                working_dataframe_end   = working_dataframe.iloc[end_point:data_length]
                working_dataframe = working_dataframe_start.append(working_dataframe_end)
            elif self.data_loss_position == 'start':
                start_point = int(np.ceil(data_length * (1-self.data_lost)))
                end_point = data_length
                working_dataframe = working_dataframe.iloc[start_point:end_point]
            elif self.data_loss_position == 'random':
                data_points_lost = int(np.floor(data_length * self.data_lost))
                keeping_samples = np.hstack((
                                        np.zeros(data_points_lost, dtype=np.bool),
                                        np.ones(data_length - data_points_lost,dtype=np.bool)
                                        ))
                self.rng.shuffle(keeping_samples)
                print(keeping_samples)
                working_dataframe = working_dataframe.iloc[np.where(keeping_samples)[0]]
            else:
                print('{} data loss method not implemented yet, keeping all data'.format(self.data_loss_position))
                start_point = 0
                end_point = data_length
            
            hourly_dataframe_out = hourly_dataframe_out.append(working_dataframe)
        
        return hourly_dataframe_out

    def sort_datasets(self, met_extracted_data, req_sites_list, var_string):
        """
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
        """
        # AG: Trim date index to be only those available in dataset: Or memory overloads and kills process.
        # Todo Doug: check OK.
        start_date = max(met_extracted_data[self._timestamp_string].min(), self.start)
        end_date = min(met_extracted_data[self._timestamp_string].max(), self.end)

        print('Start date: {}'.format(datetime.strftime(start_date, MetPostProcessor.INPUT_DATE_FORMAT)))
        print('End date: {}'.format(datetime.strftime(end_date, MetPostProcessor.INPUT_DATE_FORMAT)))
        if self.verbose > 1: print('Using date range in sort_datasets: {} to {}'.format(
            datetime.strftime(start_date, MetPostProcessor.INPUT_DATE_FORMAT),
            datetime.strftime(end_date, MetPostProcessor.INPUT_DATE_FORMAT)
        ))

        date_index = pd.date_range(start=start_date, end=end_date, freq='1H', name=self._timestamp_string)

        # add the Date index
        indexed_orig_data = met_extracted_data.set_index(self._timestamp_string)

        # define initial column for dataframe
        dataframe_columns = {var_string: np.nan}

        # empty dataframe for storing data
        full_data_out = pd.DataFrame()

        for site in req_sites_list:

            print('extract site {}'.format(site))

            work_data = indexed_orig_data[indexed_orig_data[self._site_string]==site]

            ts = pd.DataFrame(dataframe_columns, index=date_index)

            ts[var_string] = work_data[var_string]

            # as we didn't impute anything, add zero value flags
            ts['{}_flag'.format(var_string)] = 0
            ts.loc[ts[var_string].isna(),'{}_flag'.format(var_string)] = 1

            # add the site ID, and reindex
            ts[self._site_string] = site
            ts = ts.reset_index().set_index(self._columns_base)

            # copy data to new array
            full_data_out = full_data_out.append(ts)

        return full_data_out

    def organise_data(self, met_extracted_data, req_sites_temp, req_sites_pres, req_sites_dewpoint):
        """
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
        """
        print('sorting temperature data')
        met_data_out_temp = self.sort_datasets(met_extracted_data, req_sites_temp, 'temperature')
        print('sorting pressure data')
        met_data_out_pressure = self.sort_datasets(met_extracted_data, req_sites_pres, 'pressure')
        print('sorting dew point temperature')
        met_data_out_dewpoint = self.sort_datasets(met_extracted_data, req_sites_dewpoint, 'dewpoint')

        return met_data_out_temp, met_data_out_pressure, met_data_out_dewpoint


    def combine_plotting_dataset(self, met_data_in, met_data_temp, met_data_pres, met_data_rh):
        """
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
                    'temperature'         (float): hourly data
                    'temperature_flag'    (float): flag to indicate fraction of imputed data 
                                                   (1 = fully imputed, 0 = no imputed values were used)
                    'relativehumidity'         (float): hourly data
                    'relativehumidity_flag'    (float): flag to indicate fraction of imputed data 
                                                        (1 = fully imputed, 0 = no imputed values were used)
                    'pressure'         (float): hourly data
                    'pressure_flag'    (float): flag to indicate fraction of imputed data 
                                                (1 = fully imputed, 0 = no imputed values were used)
        """

        

        combined_data = met_data_temp.merge(met_data_rh, how='outer', left_index=True, right_index=True)
        combined_data = combined_data.merge(met_data_pres, how='outer', left_index=True, right_index=True)

        combined_data.rename(columns={'rel_hum':'relativehumidity','rel_hum_flag':'relativehumidity_flag'}, inplace=True)

        combined_data.sort_index(level=1,inplace=True)
        combined_data.index.rename(['time_stamp', 'sensor_name'], inplace=True)

        return combined_data


    def imputation_hourly_analysis(self,hourly_imputed_dataframe,hourly_reference_dataframe,site_list_internal):
        """
            Statistical analysis of the hourly results for the imputation of AURN data
        
        Args:
            hourly_imputed_dataframe: hourly dataset, for all measurements, as pandas.Dataframe
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
                    O3_flag       (int): flag indicating imputed data (0=original,1=imputed)
                    PM10_flag     (int):
                    PM2.5_flag    (int):
                    NO2_flag      (int):
                    NOXasNO2_flag (int):
                    SO2_flag      (int):
            hourly_reference_dataframe: hourly dataset, for all measurements, as pandas.Dataframe
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
                    O3_flag       (int): flag indicating imputed data (0=original,1=imputed)
                    PM10_flag     (int): (all these flags should be zero)
                    PM2.5_flag    (int):
                    NO2_flag      (int):
                    NOXasNO2_flag (int):
                    SO2_flag      (int):
            site_list_internal: (list, strings) a single list of required sites
        
        Returns: None
        """
        
        sns.set_style('ticks')
        sns.set_context('paper')
        sns.despine()

        note="""
        Analysis of the reliability of the imputation of Meteorological datasets
        at the site {}. This is the original hourly data.
        
        The configuration used for this test is:
        start date: {}
        end date: {}
        fraction of data removed: {}
        position in timeseries for data removal: {} 
        """
        
        spc_list = ['temperature','pressure','relativehumidity']
        mul_ind = pd.MultiIndex.from_product([site_list_internal,spc_list],names=['site_id','spc'])
        col_headers = ['kendalltau_corr','spearmanr_corr','pearsonr_corr','slope','r_squared','p_value','std_err']
        
        hourly_stat_dataset = pd.DataFrame(index=mul_ind,columns=col_headers,dtype=np.float)
        
        for site in site_list_internal:
            print('working on site: {}'.format(site))
            with PdfPages(self.pdf_file_string.format(self.out_dir,site,'hourly')) as pdf_pages:
                firstPage = plt.figure(figsize=(6,6))
                firstPage.clf()
                firstPage.text(0.5,0.5,note.format(site,self.start,self.end,self.data_lost,self.data_loss_position), 
                    transform=firstPage.transFigure, size=12, ha="center")
                pdf_pages.savefig()
                plt.close()
                for spc in spc_list: #self.species_list:
                    print('stats for species: {}'.format(spc))
                
                    data_imputed   = hourly_imputed_dataframe.loc[(slice(None),site),spc]
                    data_reference = hourly_reference_dataframe.loc[(slice(None),site),spc]
                
                    flag_imputed = hourly_imputed_dataframe.loc[(slice(None),site),'{}_flag'.format(spc)]
                
                    # keep only the data which has been imputed
                    data_imputed   = data_imputed[flag_imputed==1]
                    data_reference = data_reference[flag_imputed==1]
                
                    # remove datapoints which were NaN in the original data
                    data_imputed   = data_imputed[data_reference.notna()] 
                    data_reference = data_reference[data_reference.notna()] 
                
                    # plot scatter
                    data_combined = pd.DataFrame()
                    data_combined[spc] = data_reference
                    data_combined['{} (imputed)'.format(spc)] = data_imputed
                    
                    min_val = min(min(data_imputed),min(data_reference))
                    max_val = max(max(data_imputed),max(data_reference))
                    range_val = max_val - min_val
                    
                    k_corr, k_pval = stats.kendalltau(data_reference,data_imputed)
                    s_corr, s_pval = stats.spearmanr(data_reference,data_imputed)
                    p_corr, p_pval = stats.pearsonr(data_reference,data_imputed)
                    slope, intercept, r_value, p_value, std_err = stats.linregress(data_reference,data_imputed)
                    hourly_stat_dataset.loc[(site,spc),col_headers] = [k_corr,s_corr,p_corr,slope,r_value**2,p_value,std_err]
                    
                    
                    sns_plot = sns.jointplot(data=data_combined,x=spc,y='{} (imputed)'.format(spc),kind="reg")
                    sns_plot.ax_joint.plot(data_combined['{} (imputed)'.format(spc)],data_combined['{} (imputed)'.format(spc)], 'r-', linewidth=1)
                    sns_plot.ax_joint.set_xlim(min_val-range_val*0.05,max_val+range_val*0.05)
                    sns_plot.ax_joint.set_ylim(min_val-range_val*0.05,max_val+range_val*0.05)
                    sns_plot.ax_joint.text(min_val+range_val*0.1,max_val-range_val*0.1,'KendallTau; corr = {0:.2f}; p = {1:.2f}'.format(k_corr,k_pval))
                    
                    pdf_pages.savefig(sns_plot.fig)
                    plt.close()
                    
                    
        
        hourly_stat_dataset.to_csv(self.csv_file_string.format(self.out_dir,'hourly'), index=True, header=True, float_format=self.float_format)
        


    def imputation_daily_analysis(self,daily_imputed_dataframe,daily_reference_dataframe,site_list_internal):
        """
            Statistical analysis of the daily results for the imputation of AURN data
        
        Args:
            daily_imputed_dataframe: daily dataset, for all measurements, as pandas.Dataframe
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
            daily_reference_dataframe: daily dataset, for all measurements, as pandas.Dataframe, 
                          same layout as 'daily_imputed_dataframe'
            site_list_internal: (list, strings) a single list of required sites
        
        Returns: None
        """

        sns.set_style('ticks')
        sns.set_context('paper')
        sns.despine()

        note="""
        Analysis of the reliability of the imputation of Meteorological datasets
        at the site {}. This is the daily mean and maximum data.
        
        The configuration used for this test is:
        start date: {}
        end date: {}
        fraction of data removed: {}
        position in timeseries for data removal: {} 
        """
        spc_list = ['temperature','pressure','relativehumidity']
        stat_list = ['mean','max']
        mul_ind = pd.MultiIndex.from_product([site_list_internal,spc_list,stat_list],names=['site_id','spc','stat'])
        col_headers = ['kendalltau_corr','spearmanr_corr','pearsonr_corr','slope','r_squared','p_value','std_err']
        
        daily_stat_dataset = pd.DataFrame(index=mul_ind,columns=col_headers,dtype=np.float)

        fill_ranges = [0,0.125,0.25,0.5,1.0]
        length_ranges = len(fill_ranges)

        
        for site in site_list_internal:
            site_string = "{} [WEATHER]".format(site)
            print('working on site: {}'.format(site))
            with PdfPages(self.pdf_file_string.format(self.out_dir,site,'daily')) as pdf_pages:
                firstPage = plt.figure(figsize=(6,6))
                firstPage.clf()
                firstPage.text(0.5,0.5,note.format(site,self.start,self.end,self.data_lost,self.data_loss_position), 
                    transform=firstPage.transFigure, size=12, ha="center")
                pdf_pages.savefig()
                plt.close()
                for spc in ['temperature','pressure','relativehumidity']: #self.species_list:
                    print('stats for species: {}'.format(spc))
                
                    for stat in ['max','mean']:
                        data_imputed   = daily_imputed_dataframe.loc[(slice(None),site_string),'{}_{}'.format(spc,stat)]
                        data_reference = daily_reference_dataframe.loc[(slice(None),site_string),'{}_{}'.format(spc,stat)]
                
                        flag_imputed = daily_imputed_dataframe.loc[(slice(None),site_string),'{}_flag'.format(spc)]
                
                        flag_already_missing = daily_reference_dataframe.loc[(slice(None),site_string),'{}_flag'.format(spc)]
                        flag_plot = np.floor(1-flag_already_missing) # 0 = some missing data originally, 1 = no original missing data 
                        
                        flag_imputed = flag_imputed * flag_plot # remove all days which were originally missing some data

                        # keep the days which contain imputed data, for stat calculations
                        data_imputed_stats   = data_imputed[flag_imputed>0]
                        data_reference_stats = data_reference[flag_imputed>0]
                

                        min_val = min(min(data_imputed_stats),min(data_reference_stats))
                        max_val = max(max(data_imputed_stats),max(data_reference_stats))
                        range_val = max_val - min_val
                
                        k_corr, k_pval = stats.kendalltau(data_reference_stats,data_imputed_stats)
                        s_corr, s_pval = stats.spearmanr(data_reference_stats,data_imputed_stats)
                        p_corr, p_pval = stats.pearsonr(data_reference_stats,data_imputed_stats)
                        slope, intercept, r_value, p_value, std_err = stats.linregress(data_reference_stats,data_imputed_stats)
                        daily_stat_dataset.loc[(site,spc,stat),col_headers] = [k_corr,s_corr,p_corr,slope,r_value**2,p_value,std_err]
                
                        impute_string = '{}_{} (imputed)'.format(spc,stat)
                        data_combined = pd.DataFrame()
                        data_combined['{}_{}'.format(spc,stat)] = data_reference_stats
                        data_combined[impute_string] = data_imputed_stats
                
                        sns_plot = sns.jointplot(data=data_combined,x='{}_{}'.format(spc,stat),y=impute_string, kind='reg')
                        sns_plot.ax_joint.plot(data_combined[impute_string],data_combined[impute_string], 'r-', linewidth=1)
                        sns_plot.ax_joint.set_xlim(min_val-range_val*0.05,max_val+range_val*0.05)
                        sns_plot.ax_joint.set_ylim(min_val-range_val*0.05,max_val+range_val*0.05)
                        sns_plot.ax_joint.text(min_val+range_val*0.1,max_val-range_val*0.1,'KendallTau; corr = {0:.2f}; p = {1:.2f}'.format(k_corr,k_pval))
                        pdf_pages.savefig(sns_plot.fig)
                        plt.close()

                        for ind in range(length_ranges-1):
                        
                            # keep only the data which has been imputed
                            data_imputed_internal   = data_imputed[flag_imputed>fill_ranges[ind]]
                            data_reference_internal = data_reference[flag_imputed>fill_ranges[ind]]
                            flag_imputed_internal = flag_imputed[flag_imputed>fill_ranges[ind]]
                            
                            data_imputed_internal   = data_imputed_internal[flag_imputed_internal<=fill_ranges[ind+1]]
                            data_reference_internal = data_reference_internal[flag_imputed_internal<=fill_ranges[ind+1]]
                            
                            if not data_imputed_internal.empty:
                                # plot scatter
                                impute_string = '{}_{} (imputed, range {}-{})'.format(spc,stat,fill_ranges[ind],fill_ranges[ind+1])
                                data_combined = pd.DataFrame()
                                data_combined['{}_{}'.format(spc,stat)] = data_reference_internal
                                data_combined[impute_string] = data_imputed_internal
                
                                sns_plot = sns.jointplot(data=data_combined,x='{}_{}'.format(spc,stat),y=impute_string, kind='reg')
                                sns_plot.ax_joint.plot(data_combined[impute_string],data_combined[impute_string], 'r-', linewidth=1)
                                sns_plot.ax_joint.set_xlim(min_val-range_val*0.05,max_val+range_val*0.05)
                                sns_plot.ax_joint.set_ylim(min_val-range_val*0.05,max_val+range_val*0.05)
                                sns_plot.ax_joint.text(min_val+range_val*0.1,max_val-range_val*0.1,'KendallTau; corr = {0:.2f}; p = {1:.2f}'.format(k_corr,k_pval))
                                pdf_pages.savefig(sns_plot.fig)
                                plt.close()

                        

        daily_stat_dataset.to_csv(self.csv_file_string.format(self.out_dir,'daily'), index=True, header=True, float_format=self.float_format)
