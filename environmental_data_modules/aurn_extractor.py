try:
    import os
    import wget
    import pyreadr
    import datetime
    import pandas as pd
    import numpy as np
except:
    pass

from environmental_data_modules import EnvironmentModule, AurnModule


class AurnExtractor(EnvironmentModule, AurnModule):
    # Define 'absolute' constants
    SPECIES_LIST = ['O3', 'PM10', 'PM2.5', 'NO2', 'NOXasNO2', 'SO2']
    AVAILABLE_YEARS = [2016, 2017, 2018, 2019]
    DATE_RANGE = None  # Date range not used in Aurn Extraction, just whole years

    # Define defaults
    DEFAULT_SAVE_TO_CSV = True
    BASE_FILE_OUT = '{}/AURN_extracted{}.csv'

    def __init__(self, metadata_filename=AurnModule.DEFAULT_METADATA_FILE, metadata_url=AurnModule.DEFAULT_METADATA_URL,
                 out_dir=EnvironmentModule.DEFAULT_OUT_DIR,
                 verbose=EnvironmentModule.DEFAULT_VERBOSE):
        super(AurnExtractor, self).__init__(out_dir, verbose)
        AurnModule.__init__(self, metadata_filename=metadata_filename, metadata_url=metadata_url)
        self._base_file_out = AurnExtractor.BASE_FILE_OUT
        self.years = AurnExtractor.AVAILABLE_YEARS
        self.date_range = AurnExtractor.DATE_RANGE

    @property
    def years(self):
        return self.__years

    @years.setter
    def years(self, years):
        try:
            years = set(list(years))
        except Exception:
            raise ValueError('Years must be a list. Input: {}'.format(years))

        error_years = set(years) - set(AurnExtractor.AVAILABLE_YEARS)
        assert len(error_years) == 0, \
            "Each year must be contained in available years: {}. Error years: {}".format(
                AurnExtractor.AVAILABLE_YEARS, str(error_years))
        self.__years = years


    def extract_data(self,
                     years=AVAILABLE_YEARS,
                     site_list=AurnModule.DEFAULT_SITE_LIST,
                     save_to_csv=DEFAULT_SAVE_TO_CSV,
                     outfile_suffix=EnvironmentModule.DEFAULT_OUT_FILE_SUFFIX):

        self._outfile_suffix = outfile_suffix
        self.file_out = self._base_file_out.format(self.out_dir, self.outfile_suffix_string)
        self.years = years
        self.site_list = site_list

        # create a dataframe with the hourly dataset for all stations
        hourly_dataframe = self.extract_site_data(save_to_csv)
        hourly_dataframe = hourly_dataframe.rename(columns={AurnModule.SITE_ID_EXTRACTED: AurnModule.SITE_ID_NEW})

        # apply some filtering of negative and zero values

        for species in AurnExtractor.SPECIES_LIST:
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
        hourly_dataframe = hourly_dataframe[list(AurnModule.EXTRACTED_FILE_COLS)]

        if self.verbose > 1: print('hourly_dataframe: \n {}'.format(hourly_dataframe))
        if save_to_csv:
            # Give index a name (for saving)
            hourly_dataframe.index.name = AurnModule.INDEX_EXTRACTED
            hourly_dataframe.to_csv(self.file_out, header=True)

        return hourly_dataframe

    def extract_site_data(self, save_to_csv):

        final_dataframe = pd.DataFrame()

        for site in self._site_list:

            # select our subset of metadata for this station
            subset_df = self.metadata['AURN_metadata'][self.metadata['AURN_metadata'].site_id == site]
            station_name = subset_df['site_name'].values[0]

            print("processing site {} ({})".format(site, station_name))

            # get the list of years of data for this station
            years_process = self.define_years_to_download(subset_df, self.years)

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

            # Todo Doug: should this be removed
            # postprocessing the data set, to get daily data
            # final_dataframe = final_dataframe.append(postprocess_data(full_hourly_dataframe,site))

            if save_to_csv is True:
                # Now save the full hourly dataframe as a .csv file
                full_hourly_dataframe.to_csv(os.path.join(self.out_dir, site + '.csv'), index=False, header=True)

        return final_dataframe

        # functions for downloading data

    def define_years_to_download(self, subset_df, years):

        # Check to see if your requested years will work and if not, change it to do this.
        # Create two new columns of datetimes for earliest and latest
        datetime_start = pd.to_datetime(subset_df['start_date'].values, format='%Y/%m/%d').year
        # Problem with the end date is it could be ongoing. In which case, convert that entry into a date and to_datetime
        now = datetime.datetime.now()
        datetime_end_temp = subset_df['end_date'].values
        step = 0
        for i in datetime_end_temp:
            if i == 'ongoing':
                datetime_end_temp[step] = str(now.year) + '-' + str(now.month) + '-' + str(now.day)
            step += 1
        datetime_end = pd.to_datetime(datetime_end_temp).year

        earliest_year = np.min(datetime_start)
        latest_year = np.max(datetime_end)

        # now create list of years to process
        years_process = []

        for year in years:
            if (year <= latest_year and year >= earliest_year):
                years_process.append(year)

        return years_process

    def download_and_open_datafiles(self, site, station_name, years):

        downloaded_site_data = []

        for year in years:
            try:
                downloaded_file = site + "_" + str(year) + ".RData"
                download_url = "https://uk-air.defra.gov.uk/openair/R_data/" + downloaded_file
                print("\tdownloading file {}".format(download_url))

                # Check to see if file exists or not. Special case for current year as updates on hourly basis
                filename_path = os.path.join(self.out_dir, downloaded_file)
                if os.path.exists(filename_path):
                    print("\t\tData file already exists, will use this")
                else:
                    print("\t\tDownloading data file for ", station_name, " in ", str(year))
                    wget.download(download_url, out=str(self.out_dir))

                # Read the RData file into a Pandas dataframe
                downloaded_data = pyreadr.read_r(str(filename_path))
                # Drop non-required fields

                # Append to dataframe list
                downloaded_site_data.append(downloaded_data[site + "_" + str(year)])
            except Exception as err:
                print("\t\tCouldn't download and extract data from {} for {}. {}".format(year, station_name, err))

        return downloaded_site_data

    def load_sort_data(self, downloaded_site_data):

        final_dataframe = pd.concat(downloaded_site_data, axis=0, ignore_index=True)

        final_dataframe['datetime'] = pd.to_datetime(final_dataframe['date'])
        final_dataframe = final_dataframe.sort_values(by='datetime', ascending=True)
        # final_dataframe=final_dataframe.set_index('datetime')

        return final_dataframe

    def tidy_hourly_data(self, hourly_dataframe, site):

        columns_of_interest = ['datetime', 'O3', 'NO2', 'SO2', 'NOXasNO2', 'PM2.5', 'PM10']
        ds_columns = hourly_dataframe.columns

        # retain the data we are interested in (as not all datasets have all variables)
        columns_to_retain = set(columns_of_interest) & set(ds_columns)
        working_dataframe = hourly_dataframe[columns_to_retain]

        working_dataframe = working_dataframe.rename(columns={'datetime': 'Date'})

        # add the site as a new column, and set as part of multiindex with the date
        # site_name = "{} [AQ]".format(site)
        site_name = "{}".format(site)

        working_dataframe[AurnModule.SITE_ID_EXTRACTED] = site_name

        return working_dataframe
