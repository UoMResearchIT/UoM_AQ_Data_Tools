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

from environmental_data_modules import AurnPostProcessor

class AurnImputationTest(AurnPostProcessor):
    """
        Class for testing the imputation of the data that is extracted from the AURN server.
    """
    # testing default values
    DEFAULT_DATA_LOST = 0.5
    DEFAULT_DATA_LOSS_POSITION = 'end'
    
    
    def __init__(self, metadata_filename=AurnPostProcessor.DEFAULT_METADATA_FILE, metadata_url=AurnPostProcessor.DEFAULT_METADATA_URL,
                 out_dir=AurnPostProcessor.DEFAULT_OUT_DIR, verbose=AurnPostProcessor.DEFAULT_VERBOSE):
        """ Initialise instance of the AurnImputationTest class.
            Initialises the private class variables

            Args:
                metadata_filename: filename of the metadata used in Aurn data extraction
                metadata_url: alternative source of AURN metadata, if metadata_filename is None
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of AurnPostProcessor

        """
        super(AurnImputationTest, self).__init__(metadata_filename,metadata_url,out_dir,verbose)

        self._data_lost = AurnImputationTest.DEFAULT_DATA_LOST
        self._data_loss_position = AurnImputationTest.DEFAULT_DATA_LOSS_POSITION

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

    def imputation_test(self, in_file, date_range=None,
                site_list=AurnPostProcessor.DEFAULT_SITE_LIST,
                emep_filename=AurnPostProcessor.DEFAULT_EMEP_FILENAME,
                min_years_reference=AurnPostProcessor.DEFAULT_MIN_YEARS_REFERENCE,
                min_years=AurnPostProcessor.DEFAULT_MIN_YEARS,
                data_lost=DEFAULT_DATA_LOST, 
                data_loss_position=DEFAULT_DATA_LOSS_POSITION,
                save_to_csv=AurnPostProcessor.DEFAULT_SAVE_TO_CSV,
                outfile_suffix=''):

        """ Testing the imputation methods used for filling in missing data. Replicates the
            methods used in 'process' function, but for the defined sites will remove a
            given amount of data from the dataframe before imputation, then compare the
            imputed data with the original data, to determine accuracy of the method.
            
            The stations to be tested need to meet the requirements for 'reference' sites.
            
            Args:
                in_file:                (str) The file spec of the input file (required)
                date_range:             (list of 2 datetime) The date range of interest
                site_list:              (list of string/number) Site IDs of interest
                emep_filename:          (str) The file spec of the EMEP file to be used to help calculate #Todo Doug
                min_years_reference:    (float) The minimum number of years of data for any site that  we are going to
                                            use as a reference site later. (this cannot be less than min_years)
                min_years:              (float) The minimum number of years of data that a site must have
                data_lost:              (float) The fraction of each dataset to remove
                data_loss_position:     (str) where to lose the data (start,middle,end,random)
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
            self.date_range = [datetime.strptime(date_range[0], AurnPostProcessor.INPUT_DATE_FORMAT),
                               datetime.strptime(date_range[1], AurnPostProcessor.INPUT_DATE_FORMAT)]
        else:
            self.date_range = [self.get_available_start(), self.get_available_end()]

        self.file_out = AurnPostProcessor.BASE_FILE_OUT.format(self.out_dir, outfile_suffix)
        self._emep_data = self.load_emep_data(emep_filename)
        self.min_years = min_years
        self.min_years_reference = min_years_reference
        self.site_list = site_list
        self.data_lost = data_lost
        self.data_loss_position = data_loss_position
        self.station_data = self.metadata['AURN_metadata'][['site_id', 'latitude', 'longitude', 'site_name']]
        if self.verbose > 1: print('Station data: \n {}'.format(self.station_data))


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

        # data preparation, to create the test dataset for imputation
        hourly_test_dataframe, reference_sites = self.data_preparation(hourly_dataframe_filtered,
            site_list_internal,reference_sites)

        print('imputation of data, returning hourly data')
        hourly_imputed_dataframe = self.organise_data_imputation(
            hourly_test_dataframe, reference_sites, required_sites, site_list_internal)

        print('sorting data (no imputation), returning hourly data')
        hourly_reference_dataframe = self.organise_data(hourly_dataframe_filtered, site_list_internal)

        # calculate the daily max and mean for each station
        daily_reference_dataframe = self.combine_and_organise_mean_max(hourly_reference_dataframe)
        daily_imputed_dataframe = self.combine_and_organise_mean_max(hourly_imputed_dataframe)

        # calculate the stats for the hourly and daily data, printing out graphs of these
        self.imputation_hourly_analysis(hourly_imputed_dataframe,hourly_reference_dataframe,site_list_internal)
        self.imputation_daily_analysis(daily_imputed_dataframe,daily_reference_dataframe,site_list_internal)

        #if save_to_csv:
        #    # write this dataset to file
        #    daily_dataframe.to_csv(self.file_out, index=True, header=True, float_format='%.2f')

        #return daily_dataframe

    def data_preparation(self, hourly_dataframe, site_list_internal, reference_sites):
        """
        Prepare test data for imputation, by removing the specified amount of data from the test sites.
        
        Args:
            hourly_dataframe: hourly dataset, for all measurements, as pandas.Dataframe
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
            
            reference_sites: (dict, keys are species):
                            items: (list of strings) the siteID's for our reference sites for each `spc` 
            site_list_internal: (list, strings) a single list of required sites
        
        Returns:
            hourly_dataframe: hourly dataset, for all measurements, as pandas.Dataframe, with the required data removed
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
            reference_sites: (dict, keys are species):
                            items: (list of strings) the siteID's for our reference sites for each `spc`
                                      These will have our test sites removed, to avoid issues later 
        
        """
        
        
        return hourly_dataframe, reference_sites

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
            site_list_internal: (list, strings) a single list of required sites
        
        Returns: None
        """
        
        print('...analysis of hourly imputed data is work in progress...')

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
        
        print('...analysis of daily mean/max imputed data is work in progress...')












