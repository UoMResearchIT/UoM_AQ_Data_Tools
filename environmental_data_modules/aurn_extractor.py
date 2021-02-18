try:
    import os
    import wget
    import pyreadr
    from datetime import datetime
    import pandas as pd
    import numpy as np
    from urllib.parse import urljoin
except:
    pass

from environmental_data_modules import Extractor, AurnModule, DateYearsProcessor


class AurnExtractor(Extractor, AurnModule, DateYearsProcessor):
    """
        Class used for extracting data from the AURN server.
    """

    # Define defaults
    DEFAULT_SAVE_TO_CSV = True
    BASE_FILE_OUT = '{}/AURN_extracted{}.csv'

    def __init__(self, metadata_filename=None, metadata_url=AurnModule.DEFAULT_METADATA_URL,
                 out_dir=Extractor.DEFAULT_OUT_DIR,
                 verbose=Extractor.DEFAULT_VERBOSE):
        """ Initialise instance of the AurnExtractor class.
            Initialises the private class variables

            Args:
                metadata_filename: filename of the metadata used in Aurn data extraction
                metadata_url: alternative source of AURN metadata, if metadata_filename is None
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of AurnExtractor

        """
        super(AurnExtractor, self).__init__(out_dir, verbose)
        AurnModule.__init__(self, metadata_filename=metadata_filename, metadata_url=metadata_url)
        DateYearsProcessor.__init__(self)
        self._base_file_out = AurnExtractor.BASE_FILE_OUT


    def extract_data(self,
                     years=DateYearsProcessor.get_available_years(),
                     site_list=AurnModule.DEFAULT_SITE_LIST,
                     save_to_csv=DEFAULT_SAVE_TO_CSV,
                     outfile_suffix=Extractor.DEFAULT_OUT_FILE_SUFFIX,
                     species_list=AurnModule.SPECIES_LIST_EXTRACTED):
        """ Extracts the AURN data for the given years and sites from the Rdata files
            downloaded from the AURN server.

            Args:
                years:              (list of integer) The years of interest
                site_list:          (list of numbers/strings) The site IDs of interest
                save_to_csv:        (boolean) Whether to save the output dateframes to CSV file(s)
                outfile_suffix:     (string) The suffix to appended to the end of output file names.
                species_list:       (list of string).

            Returns:
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
        """

        self._outfile_suffix = outfile_suffix
        self.file_out = self._base_file_out.format(self.out_dir, self.outfile_suffix_string)
        self.years = years
        self.site_list = site_list
        self.species_list = species_list

        # create a dataframe with the hourly dataset for all stations
        hourly_dataframe = self.extract_site_data(save_to_csv)

        # apply some filtering of negative and zero values

        for species in self.species_list:
            print('filtering {}:'.format(species))
            try:
                print('\t{} has {} positive values'.format(species,
                                                           len(hourly_dataframe.loc[hourly_dataframe[species] > 0.0])))
                print('\t{} has {} NaNs'.format(species, len(hourly_dataframe.loc[hourly_dataframe[species].isna()])))
                print('\t{} has {} negative or zero values that will be replaced with NaNs'.format(species, len(
                    hourly_dataframe.loc[hourly_dataframe[species] <= 0.0])))
                hourly_dataframe.loc[hourly_dataframe[species] <= 0.0, species] = np.nan
            except:
                print('\t{} has  no values'.format(species))

        # Put column names in correct order
        hourly_dataframe = hourly_dataframe[list(AurnModule.NEW_FILE_COLS)]

        if self.verbose > 1: print('hourly_dataframe: \n {}'.format(hourly_dataframe))
        if save_to_csv:
            # Give index a name (for saving)
            hourly_dataframe.index.name = AurnModule.INDEX_EXTRACTED
            hourly_dataframe.to_csv(self.file_out, header=True)

        return hourly_dataframe

    def extract_site_data(self, save_to_csv):

        """ Organises the downloading and extraction of the AURN dataset.

            Returns:
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
        """

        final_dataframe = pd.DataFrame()

        for site in self._site_list:

            # select our subset of metadata for this station
            subset_df = self.metadata['AURN_metadata'][self.metadata['AURN_metadata'][AurnModule.SITE_ID_AURN_METADATA] == site]
            station_name = subset_df['site_name'].values[0]

            print("processing site {} ({})".format(site, station_name))

            # get the list of years of data for this station
            years_process = self.define_years_to_download(subset_df)

            # if the list of years is empty, then skip this site
            if not years_process:
                print("\t\tNo data for years of interest, skipping this site")
                continue

            # download the datasets
            downloaded_site_data = self.download_and_open_datafiles(site, station_name, years_process)

            # if we couldn't download the data, skip this site
            if len(downloaded_site_data) == 0:
                print("\t\tNo data could be downloaded for {}".format(station_name))
                continue

            # combine and sort data
            full_hourly_dataframe = self.load_sort_data(downloaded_site_data)

            # tidy data, adding site column and removing excess columns
            full_hourly_dataframe = self.tidy_hourly_data(full_hourly_dataframe, site)

            # return the full hourly dataset
            final_dataframe = final_dataframe.append(full_hourly_dataframe)

            if save_to_csv is True:
                # Now save the full hourly dataframe as a .csv file
                full_hourly_dataframe.to_csv(os.path.join(self.out_dir, site + '.csv'), index=False, header=True)

        return final_dataframe

    def define_years_to_download(self, subset_df):

        """ Checks what years of data are available for the selected station.
        
            Args:
                subset_df: AURN metadata for given station, as a pandas.Dataframe
                    Required Columns:
                        start_date   (str): YYYY-MM-DD, starting date for measurements
                        end_date     (str): end date, or 'ongoing'
                
            Returns:
                years_process  (list of ints): the years to process data for this site
        """

        # extract starting years for measurement site (for all pollutants measured there)
        datetime_start = pd.to_datetime(subset_df['start_date'].values, format='%Y/%m/%d').year
        # extract the end years for measurement site (for all pollutants measured there)- 
        # if this is given as 'ongoing' then use the current year
        now = datetime.now()
        datetime_end_temp = subset_df['end_date'].values
        step = 0
        for i in datetime_end_temp:
            if i == 'ongoing':
                datetime_end_temp[step] = str(now.year) + '-' + str(now.month) + '-' + str(now.day)
            step += 1
        datetime_end = pd.to_datetime(datetime_end_temp).year

        # get single values for start and end years
        earliest_year = np.min(datetime_start)
        latest_year = np.max(datetime_end)

        # compare years requested with years available, to create a list of years that
        # can actually be processed for this site
        years_process = []

        for year in self.years:
            if year <= latest_year and year >= earliest_year:
                years_process.append(year)

        return years_process

    def download_and_open_datafiles(self, site, station_name, years):

        """ Downloads and opens the AURN datafiles
        
            Args:
                site           (str): siteID, used in the file names
                station_name   (str): used for diagnostic messages
                years (list of ints): list of the years to download
                         
            Returns:
                downloaded_site_data   (list): list of the dataframes that have been downloaded
        """

        downloaded_site_data = []

        for year in years:
            try:
                downloaded_file = "{}_{}.RData".format(site, str(year))
                download_url = urljoin(AurnExtractor.DEFAULT_DOWNLOAD_RDATA_URL, downloaded_file)
                print("\tdownloading file {}".format(download_url))

                # Check to see if file exists or not. Special case for current year as updates on hourly basis
                filename_path = os.path.join(self.out_dir, downloaded_file)
                if os.path.exists(filename_path):
                    print("\t\tData file already exists, will use this")
                else:
                    print("\t\tDownloading data file for {} in {}".format(station_name, str(year)))
                    wget.download(download_url, out=str(self.out_dir))

                # Read the RData file into a Pandas dataframe
                downloaded_data = pyreadr.read_r(str(filename_path))

                # Append to dataframe list
                downloaded_site_data.append(downloaded_data[site + "_" + str(year)])
            except Exception as err:
                print("\t\tCouldn't download and extract data from {} for {}. {}".format(year, station_name, err))

        return downloaded_site_data

    def load_sort_data(self, downloaded_site_data):
        """ Sorts the AURN data, combining these into a single dataframe, and converting
            the date strings into datetime objects.
        
            Args:
                downloaded_site_data   (list): list of the dataframes that have been downloaded

            Returns:
                final_dataframe: hourly dataset, for the given site, as pandas.Dataframe
                    Index: none
                    Required Columns:
                        Date   (datetime object):
                    Optional Columns:
                        O3       (float):
                        PM10     (float):
                        PM2.5    (float):
                        NO2      (float):
                        NOXasNO2 (float):
                        SO2      (float):
                        (extra columns, for other pollutants)
        """

        final_dataframe = pd.concat(downloaded_site_data, axis=0, ignore_index=True)
        final_dataframe[self._timestamp_string] = pd.to_datetime(final_dataframe[AurnModule.DATE_EXTRACTED])
        final_dataframe.drop(columns=[AurnModule.DATE_EXTRACTED], inplace=True)
        final_dataframe = final_dataframe.sort_values(by=self._timestamp_string, ascending=True)

        return final_dataframe

    def tidy_hourly_data(self, hourly_dataframe, site_name):
        """ Removes any unneeded pollutant data, and adds a column with the siteID.
        
            Args:
                hourly_dataframe: hourly dataset, for the given site, as pandas.Dataframe
                    Index: none
                    Required Columns:
                        Date   (datetime object):
                    Optional Columns:
                        O3       (float):
                        PM10     (float):
                        PM2.5    (float):
                        NO2      (float):
                        NOXasNO2 (float):
                        SO2      (float):
                        (extra columns, for other pollutants)

            Returns:
                working_dataframe: hourly dataset, for the given site, as pandas.Dataframe
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
        columns_of_interest = [self._timestamp_string] + self.species_list
        ds_columns = hourly_dataframe.columns

        # retain the data we are interested in (as not all datasets have all variables)
        columns_to_retain = set(columns_of_interest) & set(ds_columns)
        working_dataframe = hourly_dataframe[columns_to_retain].copy()
        working_dataframe.loc[:, self._site_string] = site_name

        return working_dataframe
