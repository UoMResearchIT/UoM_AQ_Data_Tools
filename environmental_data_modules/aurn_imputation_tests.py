try:
    import os, errno
    from datetime import datetime
    import pandas as pd
    import numpy as np
    from pathlib import Path

    from sklearn.experimental import enable_iterative_imputer
    from sklearn.impute import IterativeImputer
    from sklearn.linear_model import BayesianRidge
    from sklearn import preprocessing
    
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    from scipy import stats

    sns.set_theme()
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
    DEFAULT_CHECK_SITES = False
    DEFAULT_STAT_DIR = './'
    
    BASE_IMPUTED_STATS_PDF_FILE = '{}/{}_{}_imputed_comparison.pdf'
    BASE_IMPUTED_STATS_CSV_FILE = '{}/aurn_{}_correlation_stats.csv'
    DEFAULT_FLOAT_FORMAT = '%.4f'
    
    def __init__(self, metadata_filename=AurnPostProcessor.DEFAULT_METADATA_FILE, metadata_url=AurnPostProcessor.DEFAULT_METADATA_URL,
                 out_dir=AurnPostProcessor.DEFAULT_OUT_DIR, verbose=AurnPostProcessor.DEFAULT_VERBOSE, stat_dir=AurnImputationTest.DEFAULT_STAT_DIR):
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
        self.check_sites = AurnImputationTest.DEFAULT_CHECK_SITES
        self.rng = np.random.RandomState(0)
        self.stat_dir = stat_dir

    @property
    def stat_dir(self):
        return self.__stat_dir

    @stat_dir.setter
    def stat_dir(self, dir_name):
        try:
            dir_name = str(dir_name)
        except ValueError as err:
            raise err
        try:
            os.makedirs(dir_name)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise ValueError("Directory name {} cannot be created.".format(dir_name))
        self.__stat_dir = dir_name

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
                outfile_suffix='',check_sites=DEFAULT_CHECK_SITES,
                species_list=AurnPostProcessor.SPECIES_LIST_EXTRACTED):

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
                check_sites:            (boolean) If True then routine will list appropriate stations
                                            to use for the imputation tests, then exit
                species_list:            (list of strings) list of chemical species to test

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
        self.species_list = species_list
        self.check_sites = check_sites
        self.station_data = self.metadata['AURN_metadata'][['site_id', 'latitude', 'longitude', 'site_name']]
        if self.verbose > 1: print('Station data: \n {}'.format(self.station_data))
        self.pdf_file_string = AurnImputationTest.BASE_IMPUTED_STATS_PDF_FILE
        self.csv_file_string = AurnImputationTest.BASE_IMPUTED_STATS_CSV_FILE
        self.float_format = AurnImputationTest.DEFAULT_FLOAT_FORMAT


        # load and prepare the hourly dataset
        hourly_dataframe = self.load_aurn_data(in_file)

        print('filter for minimum data lengths, and reduce dataset to only stations of interest')
        hourly_dataframe_filtered, reference_sites, required_sites, site_list_internal = \
            self.site_list_and_preparation(hourly_dataframe)
        
        if len(hourly_dataframe_filtered.index) == 0:
            print('Exiting post-processing: Metadata is empty after initial filtering processes')
            return

        print('data preparation, to create the test dataset for imputation')
        hourly_test_dataframe = self.data_preparation(hourly_dataframe_filtered, 
                                               reference_sites, site_list_internal)

#        return

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

    def site_list_and_preparation(self,hourly_dataframe):
        """
        Wrapper for the list_required_and_reference_sites routine. This will list sites 
        based on the given minimum and reference year requirements, then determine which
        sites for all species of interest fit the reference year requirements. If none do,
        or if 'check_sites' flag is True, then potential sites of use will be listed, and 
        the program exited.
        
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
            
        Returns:
            hourly_dataframe_filtered: pandas dataframe, as above, containing hourly dataset for only 
                               the reference station datasets
            reference_sites_out: (dict, keys are species):
                            items: (list of strings) the siteID's for our reference sites for each `spc`
                                                     minus the sites in the `site_working` list
            required_sites: (dict, keys are species):
                            items: (list of strings) required sites for `spc` (not used later?)
            site_working_list: (list, strings) a single list of required sites for the imputation tests
        """

        df_part_filtered, reference_sites, required_sites, site_list_internal = \
            self.list_required_and_reference_sites(hourly_dataframe)

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
            if len(reference_sites_out[spc]) == 0:
                print('there are no reference sites for species {}, please set a reduced site list'.format(spc))
                return pd.DataFrame(), [], [], []

        # success! return the filtered dataframe, and our lists of sites
        return hourly_dataframe_filtered, reference_sites_out, required_sites, site_working_list

    def data_preparation(self, hourly_dataframe, reference_sites, site_list_internal):
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
                                                     minus the sites in the `site_working` list
            site_list_internal: (list, strings) a single list of required sites
        
        Returns:
            hourly_dataframe_out: hourly dataset, for all measurements, as pandas.Dataframe, with the required data removed
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

        date_index = pd.date_range(start=self.start, end=self.end, freq='1H', name=self._timestamp_string)
        output_dataframe = pd.DataFrame()

        hourly_dataframe_internal = hourly_dataframe_filtered.set_index(self._timestamp_string)
        spc_list = self.species_list

        if self.verbose > 1: print('1. Site list internal: ', site_list_internal)
        for site in site_list_internal:
            if self.verbose > 1: print('2. Site: ', site)

            # create new dataframe, with the dates that we are interested in
            working_hourly_dataframe = pd.DataFrame([], index=date_index)
            working_hourly_dataframe[self._site_string] = site

            # copy these to a new dataframe
            working_hourly_dataframe[spc_list] = \
                hourly_dataframe_internal[hourly_dataframe_internal[self._site_string] == site][spc_list]

            # copy imputed data of interest into copy of original dataframe (without EMEP and neighbouring sites)
            for spc in spc_list:
                working_hourly_dataframe['{}_flag'.format(spc)] = 0
                working_hourly_dataframe.loc[working_hourly_dataframe[spc].isna(),'{}_flag'.format(spc)] = 1

            # append data to the output dataframe
            output_dataframe = output_dataframe.append(working_hourly_dataframe)

        output_dataframe = output_dataframe.reset_index().set_index([self._timestamp_string,self._site_string])
        return(output_dataframe)


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
        Analysis of the reliability of the imputation of AURN datasets
        at the site {}. This is the original hourly data.
        
        The configuration used for this test is:
        start date: {}
        end date: {}
        fraction of data removed: {}
        position in timeseries for data removal: {} 
        """

        mul_ind = pd.MultiIndex.from_product([site_list_internal,self.species_list],names=[self._site_string,'spc'])
        col_headers = ['kendalltau_corr','spearmanr_corr','pearsonr_corr','slope','r_squared','p_value','std_err']
        
        hourly_stat_dataset = pd.DataFrame(index=mul_ind,columns=col_headers,dtype=np.float)
        
        for site in site_list_internal:
            print('working on site: {}'.format(site))
            with PdfPages(self.pdf_file_string.format(self.stat_dir,site,'hourly')) as pdf_pages:
                firstPage = plt.figure(figsize=(6,6))
                firstPage.clf()
                firstPage.text(0.5,0.5,note.format(site,self.start,self.end,self.data_lost,self.data_loss_position), 
                    transform=firstPage.transFigure, size=12, ha="center")
                pdf_pages.savefig()
                plt.close()
                for spc in self.species_list:
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
                    
        
        hourly_stat_dataset.to_csv(self.csv_file_string.format(self.stat_dir,'hourly'), index=True, header=True, float_format=self.float_format)


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
        Analysis of the reliability of the imputation of AURN datasets
        at the site {}. This is the daily mean and maximum data.
        
        The configuration used for this test is:
        start date: {}
        end date: {}
        fraction of data removed: {}
        position in timeseries for data removal: {} 
        """
        
        stat_list = ['mean','max']
        mul_ind = pd.MultiIndex.from_product([site_list_internal,self.species_list,stat_list],names=[self._site_string,'spc','stat'])
        col_headers = ['kendalltau_corr','spearmanr_corr','pearsonr_corr','slope','r_squared','p_value','std_err']

        daily_stat_dataset = pd.DataFrame(index=mul_ind,columns=col_headers,dtype=np.float)
        
        fill_ranges = [0,0.125,0.25,0.5,1.0]
        length_ranges = len(fill_ranges)
        
        for site in site_list_internal:
            site_string = "{} [AQ]".format(site)
            print('working on site: {}'.format(site))
            with PdfPages(self.pdf_file_string.format(self.stat_dir,site,'daily')) as pdf_pages:
                firstPage = plt.figure(figsize=(6,6))
                firstPage.clf()
                firstPage.text(0.5,0.5,note.format(site,self.start,self.end,self.data_lost,self.data_loss_position), 
                    transform=firstPage.transFigure, size=12, ha="center")
                pdf_pages.savefig()
                plt.close()
                for spc in self.species_list:
                    print('stats for species: {}'.format(spc))
                
                    for stat in stat_list:
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
                
                        sns_plot = sns.jointplot(data=data_combined,x='{}_{}'.format(spc,stat),y=impute_string, kind="reg")
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
                
                                sns_plot = sns.jointplot(data=data_combined,x='{}_{}'.format(spc,stat),y=impute_string,kind="reg")
                                sns_plot.ax_joint.plot(data_combined[impute_string],data_combined[impute_string], 'r-', linewidth=1)
                                sns_plot.ax_joint.set_xlim(min_val-range_val*0.05,max_val+range_val*0.05)
                                sns_plot.ax_joint.set_ylim(min_val-range_val*0.05,max_val+range_val*0.05)
                                sns_plot.ax_joint.text(min_val+range_val*0.1,max_val-range_val*0.1,'KendallTau; corr = {0:.2f}; p = {1:.2f}'.format(k_corr,k_pval))
                                pdf_pages.savefig(sns_plot.fig)
                                plt.close()


        daily_stat_dataset.to_csv(self.csv_file_string.format(self.stat_dir,'daily'), index=True, header=True, float_format=self.float_format)












