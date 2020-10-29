#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
import metpy.calc as mpcalc
from metpy.units import units
import geopy.distance as distance
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn import preprocessing

from environmental_data_workflows.environment_workflow import EnvironmentWorkflow

class PostProcessor(EnvironmentWorkflow):
    __metaclass__ = ABCMeta

    DEFAULT_IMPUTE_DATA = True
    DEFAULT_PRINT_STATS = True
    DEFAULT_FILE_IN = 'data_met/temp_rh_press_dewtemp_2016-2019.csv'
    DEFAULT_FILE_OUT = 'daily_mean_max_temp_RH_pres.csv'


    def __init__(self, out_dir=DEFAULT_OUT_DIR, verbose=DEFAULT_VERBOSE):
        super(PostProcessor, self).__init__(out_dir, verbose)

        self.impute_data = MetProcessor.DEFAULT_IMPUTE_DATA
        self.print_stats = MetProcessor.DEFAULT_PRINT_STATS
        self.file_in = MetProcessor.DEFAULT_FILE_IN
        self.file_out = MetProcessor.DEFAULT_FILE_OUT


    #%% station geographic routines
    def calc_station_distances(self, stations_in, stat_location):

        station_distances = pd.DataFrame(index=stations_in.index)
        station_distances['Distance'] = np.nan

        for index, row in stations_in.iterrows():
            new_location = (row['Latitude'],row['Longitude'])

            station_distances.loc[index]['Distance'] = distance.distance(stat_location,new_location).km

        return(station_distances)



class MetPostProcessor(PostProcessor):

    DEFAULT_COL_LIST = ['date', 'siteID', 'temperature', 'rh', 'pressure', 'dewtemp']
    DEFAULT_STATION_DROP_LIST = [117]
    DEFAULT_MIN_TEMPERATURE = -20
    DEFAULT_STATIONS_FILENAME = "station_data/station_data_clean.csv"
    DEFAULT_REFERENCE_STATION_NUMBER = 5

    DEFAULT_IMPUTER_RANDOM_STATE = 0
    DEFAULT_IMPUTER_ADD_INDICATOR = True
    DEFAULT_IMPUTER_INITIAL_STRATEGY = 'mean'
    DEFAULT_IMPUTER_MAX_ITER = 300
    DEFAULT_IMPUTER_ESTIMATOR = BayesianRidge()

    DEFAULT_QT_OUTPUT_DISTRIBUTION = 'normal'
    DEFAULT_QT_RANDOM_STATE = 0

    def __init__(self, out_dir=DEFAULT_OUT_DIR, verbose=DEFAULT_VERBOSE):
        super(MetProcessor, self).__init__(out_dir, verbose)

        self._col_list = MetExtractor.DEFAULT_COL_LIST
        self.station_drop_list = MetProcessor.DEFAULT_STATION_DROP_LIST
        self.min_temperature = MetProcessor.DEFAULT_MIN_TEMPERATURE
        self.stations = pd.read_csv(DEFAULT_STATIONS_FILENAME)
        self.stations = self.stations.set_index('Station')
        reference_station_number = MetProcessor.DEFAULT_REFERENCE_STATION_NUMBER

    def post_process(self, file_in, date_range=EnvironmentWorkflow.DEFAULT_DATE_RANGE,
                     file_out=PostProcessor.DEFAULT_FILE_OUT, station_drop_list=DEFAULT_STATION_DROP_LIST,
                     min_temperature=DEFAULT_MIN_TEMPERATURE, reference_station_number=DEFAULT_REFERENCE_STATION_NUMBER,
                     impute_data=PostProcessor.DEFAULT_IMPUTE_DATA, print_stats=PostProcessor.DEFAULT_PRINT_STATS,
                     random_state=DEFAULT_IMPUTER_RANDOM_STATE, add_indicator=DEFAULT_IMPUTER_ADD_INDICATOR,
                     inital_strategy=DEFAULT_IMPUTER_INITIAL_STRATEGY,
                     max_iter=DEFAULT_IMPUTER_MAX_ITER, estimator=DEFAULT_IMPUTER_ESTIMATOR,
                     output_distribution=DEFAULT_QT_OUTPUT_DISTRIBUTION, random_state=DEFAULT_QT_RANDOM_STATE):
        self.imputer = IterativeImputer(random_state=random_state, add_indicator=add_indicator,
                                   initial_strategy=inital_strategy, max_iter=max_iter, verbose=self.verbose,
                                   estimator=estimator)

        # set the power transform options
        pt = preprocessing.QuantileTransformer(output_distribution, random_state=random_state)

        print('loading met data file')
        met_data = self.load_met_data(file_in)

        print('dropping duplicate values and unwanted stations')
        met_data = self.find_and_drop_duplicates_and_unwanted_stations(met_data, station_drop_list)

        print('dropping single daily measurement stations')
        met_data = self.drop_single_daily_measurement_stations(met_data, print_stats)

        print('filtering to remove unrealistically low temperatures')
        met_data = self.remove_low_temperature_data(met_data, min_temperature)

        print('filter for minimum data lengths, and reduce dataset to only stations of interest')
        met_data, reference_sites, req_sites_temp, req_sites_pres, req_sites_dewtemp = \
            self.list_required_and_reference_sites(met_data)

        if impute_data:
            print('imputation of data, returning hourly data')
            met_data_temp, met_data_pres, met_data_dewtemp = self.organise_data_imputation(met_data,
                date_range[0], date_range[1], stations, reference_sites, req_sites_temp,
                req_sites_pres, req_sites_dewtemp, reference_station_number)
        else:
            print('sorting data (no imputation), returning hourly data')
            met_data_temp, met_data_pres, met_data_dewtemp = self.organise_data(met_data, start_date, end_date, \
                                                                           req_sites_temp, req_sites_pres,
                                                                           req_sites_dewtemp)

        print('calculation of relative humidity from temperature and dew point temperature')
        met_data_rh = self.rh_calculations(met_data_temp, met_data_dewtemp, print_stats, met_data)

        # calculate the daily max and mean for each station
        met_data_hourly = self.combine_and_organise_mean_max(met_data_temp, met_data_pres, met_data_rh)

        # write data to file
        met_data_hourly.to_csv(file_out, index=True, header=True, float_format='%.2f')



    def load_met_data(self, file_in):

        print('    load data file')
        met_data = pd.read_csv(file_in,usecols=self._col_list,skiprows=1,engine='python',sep=',\s*',na_values='None')

        print('    correct date string')
        met_data['Date'] = met_data['date'].apply(calc_date)
        met_data = met_data.drop(columns='date')

        return(met_data)


    #%% function for writing out some information about the data count stats

    def print_data_count_stats(self, dc_in):

        print('total temperature daily data count is: {}'.format(dc_in.count().values[0]))
        print('# data points per day, total daily data point counts')
        print('      temperature, rel hum, pressure, dew point temp')
        for dpoint in range(0,49):
            dcount = dc_in[dc_in.temperature==dpoint].count().values[0]
            ecount = dc_in[dc_in.rh==dpoint].count().values[0]
            fcount = dc_in[dc_in.pressure==dpoint].count().values[0]
            gcount = dc_in[dc_in.dewtemp==dpoint].count().values[0]
            print('{0:4},{1:7},{2:7},{3:7},{4:7}'.format(dpoint,dcount,ecount,fcount,gcount))


    #%% functions for finding and removing unwanted data

    def find_and_drop_duplicates_and_unwanted_stations(self, met_data_in, station_list):
        # Pull out all duplicated values
        met_duplicates = met_data_in[met_data_in.duplicated(subset=['Date','siteID'],keep=False)]

        # Split these into those with, and without, pressure data (SYNOP will have pressure data)
        #   We will keep (most of) the readings with pressure data, and will drop (most of) the
        #   readings without pressure data.
        met_dup_with_pres = met_duplicates[met_duplicates['pressure'].isna()==False]
        met_dup_no_pres = met_duplicates[met_duplicates['pressure'].isna()==True]

        # Get the 2nd data points which are duplicated *and* where both have pressure data (these are only a few points)
        #   These will be the few readings with pressure data that we drop, so we set keep to 'first' so that we
        #   get the indexes for the 'last' values.
        met_dup_with_pres_dups = met_dup_with_pres[met_dup_with_pres.duplicated(subset=['Date','siteID'],keep='first')]

        # Find the duplicates where there is no pressure data for either, and save the first value (<10,000 points)
        #   In this case we also want to preserve the first data points, but as we are going to drop the
        #   indexes we extract from our list of indexes to throw away we want to set keep to 'last' in this instance.
        met_dup_no_pres_dups = met_dup_no_pres[met_dup_no_pres.duplicated(subset=['Date','siteID'],keep='last')]
        met_dup_no_pres_reduced = met_dup_no_pres.drop(index=met_dup_no_pres_dups.index)


        # Build an index of the duplicated datapoints to drop - starting with the data with no pressure readings
        #   that will be dropped, then appending the indexes of the few data points with pressure data that we don't want.
        indexes_to_drop = met_dup_with_pres_dups.index
        indexes_to_drop = indexes_to_drop.append(met_dup_no_pres_reduced.index)


        # Append to this list the indexes of all station data that we are dropping completely
        for station in station_list:
            station_drop_indexes = met_data_in[met_data_in['siteID']==station].index
            indexes_to_drop = indexes_to_drop.append(station_drop_indexes)


        # Finally drop all the data that is unwanted
        met_data_reduced = met_data_in.drop(index=indexes_to_drop)

        return(met_data_reduced)


    def drop_single_daily_measurement_stations(self, met_data_in, print_stats):
        # group the data by date, and count the readings per day
        tempgroups = met_data_in.groupby(['siteID',pd.Grouper(key='Date', freq='1D')])
        data_counts = tempgroups.count()

        # some diagnostic output, if required
        if print_stats:
            print_data_count_stats(data_counts)

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

        if print_stats:
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
        if print_stats:
            print('new stats for the reduced data:')
            tempgroups = met_data_reduced.groupby(['siteID',pd.Grouper(key='Date', freq='1D')])
            data_counts = tempgroups.count()
            print_data_count_stats(data_counts)


        return(met_data_reduced)


    #%% function for getting two lists of stations, one for required site, one for reference sites

    def station_listing(self, var_string, met_data_in, min_years=1, reference_num_years=3.5):
        '''
        arguments:
            var_string:
                label for variable of interest
            met_data_in:
                measurement data set
            min_years (default 1):
                minimum number of years of data that a site must have
            reference_num_years (default 3.5):
                minimum number of years of data for any site that we
                are going to use as a reference site later

        returns:
            required_site_list:
                list of sites with a data count > min_years
            reference_site_list:
                list of sites with a data count > useful_num_years
        '''

        site_list_interior = met_data_in['siteID'].unique()

        required_site_list=[]
        reference_site_list=[]

        for site in site_list_interior:
            metsite = met_data_in[met_data_in['siteID']==site]
            try:
                measurement_num = len(metsite[metsite[var_string].notna()])
            except:
                measurement_num = 0
            if(measurement_num > min_years*365*24):
                required_site_list.append(site)
            if(measurement_num > reference_num_years*365*24):
                reference_site_list.append(site)

        return(required_site_list,reference_site_list)


    def list_required_and_reference_sites(self, met_data_in):
        print('    get the lists of required and reference stations for each measurement variable')
        req_sites_temp,    reference_sites_temp     = station_listing('temperature',met_data_in)
        req_sites_pres,    reference_sites_pres     = station_listing('pressure',met_data_in)
        req_sites_dewtemp, reference_sites_dewtemp  = station_listing('dewtemp',met_data_in)


        # find a unified list of useful sites for all our measurements
        reference_sites = [x for x in reference_sites_dewtemp if x in reference_sites_pres]

        # get a list of all sites which are required for at least one measurement set
        required_sites = list(dict.fromkeys(req_sites_temp + req_sites_pres + req_sites_dewtemp))

        print('there are {} required sites, and {} reference sites'.format(len(required_sites),len(reference_sites)))

        met_data_filtered = met_data[met_data['siteID'].isin(required_sites)]

        return(met_data_filtered,reference_sites,req_sites_temp,req_sites_pres,req_sites_dewtemp)


    #%% functions for calculating RH from temperature and dew point temperature data

    def rh_calculations(self, met_data_temp,met_data_dewtemp,print_stats,met_data):
        # merge the two input datasets, dropping indexes which are not in both
        met_data_in = met_data_temp.merge(met_data_dewtemp,how='inner',left_index=True,right_index=True)

        # create output data frame
        met_data_out = pd.DataFrame(index=met_data_in.index)

        met_data_out['rh'] = mpcalc.relative_humidity_from_dewpoint(met_data_in['temperature'].values * units.degC, \
                                      met_data_in['dewtemp'].values * units.degC) * 100.0


        met_data_out['rh.flag'] = met_data_in[['temperature.flag','dewtemp.flag']].values.max(1)

        # plot the distribution of the calculated RH difference from measured RH
        if print_stats:
            met_data_internal = met_data.set_index(['Date','siteID'])

            met_data_internal['rh2'] = met_data_out['rh']

            met_data_internal['rhdiff'] = met_data_internal['rh2'] - met_data_internal['rh']

            ax = met_data_internal['rhdiff'].hist(bins=100)
            ax.semilogy()
            ax.set_xlabel('RH_new - RH')


        return(met_data_out)


    #%% station geographic routines

    def get_station_distances(self, stations_in,site_in,useful_sites_in):

        station_location = (stations_in.loc[site_in]['Latitude'],stations_in.loc[site_in]['Longitude'])
        station_distances = calc_station_distances(stations_in=stations_in.loc[useful_sites_in], \
                                                   stat_location=station_location)

        # sort by distance, then drop any station which is the same location as our site of interest
        station_distances = station_distances.sort_values(by='Distance',ascending=True)
        station_distances[station_distances.Distance==0]=np.nan
        station_distances = station_distances.dropna()

        return(station_distances)




    #%% functions for imputation of the datasets

    def transform_and_impute_data(self, df_in,pt,imputer):

        # define the method we wish to use
        #pt = preprocessing.PowerTransformer(method='box-cox', standardize=False)

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

        # power transformer fitting and transforming
        pt.fit(df_work.dropna())
        np_out = pt.transform(df_work)

        # impute the missing values in this new dataframe
        imputer.fit(np_out)
        imp_out = imputer.transform(np_out)

        # apply the inverse transformation for our datasets (leaving out the indicator flags)
        np_inv = pt.inverse_transform(imp_out[:,:np_out.shape[1]])

        # copy the transformed values to a new dataframe
        df_out = df_in.copy(deep=True)
        for pos,col in enumerate(col_save):
            pos_out = list(cols).index(col)
            df_out.iloc[:,pos_out] = np_inv[:,pos]

        return(df_out)



    def get_full_datasets(self, imputer,pt,orig_data_in,req_sites_list,useful_sites_list,stations,
                          date_index,station_number,
                          station_list_string,var_string):

        # add the Date index
        indexed_orig_data = orig_data_in.set_index('Date')

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

                station_distances = get_station_distances(stations,site,useful_sites_list)
                # get data for the 5 closest stations:
                for ii in range(0,station_number):
                    ts[station_list_string[ii]] = \
                        indexed_orig_data[indexed_orig_data.siteID==station_distances.index[ii]][var_string]


                # run the imputation process
                imputed_hourly_dataframe = transform_and_impute_data(ts,pt=pt,imputer=imputer)

                # drop the extra columns of station data
                ts = ts.drop(columns=station_list_string)

                # add a flag to denote the values which have been imputed
                #   then replace the original data with the imputed data
                ts['{}.flag'.format(var_string)] = ts[var_string].isna() * 1
                ts[var_string] = imputed_hourly_dataframe[var_string]

            else:
                # if we didn't impute anything, add zero value flags
                ts['{}.flag'.format(var_string)] = 0



            # now extract the daily max and mean data
            #out_data = extract_mean_max(ts,data_template,site,var_in_string,var_out_string)

            # add the site ID, and reindex
            ts['siteID'] = site
            ts = ts.reset_index().set_index(['Date','siteID'])

            # copy data to new array
            full_data_out = full_data_out.append(ts)


        return(full_data_out)



    def organise_data_imputation(self, met_data_in,start_date_in,end_date_in,stations,\
                                     reference_sites, req_sites_temp, req_sites_pres, req_sites_dewtemp, \
                                     reference_station_number):

        date_index = pd.date_range(start=start_date, end=end_date, freq='1H', name='Date')

        station_list = [ 'station{}'.format(x+1) for x in range(0,reference_station_number) ]

        print('imputing temperature data')
        met_data_out_temp = get_full_datasets(self.imputer, self.pt, met_data_in, req_sites_temp,
                                              reference_sites,stations, date_index, reference_station_number,
                                              station_list, 'temperature')

        print('imputing pressure data')
        met_data_out_pressure = get_full_datasets(self.imputer, self.pt, met_data_in, req_sites_pres, reference_sites,
                                                  stations, date_index,reference_station_number, station_list,'pressure')

        print('imputing dew point temperature')
        met_data_out_dewtemp = get_full_datasets(self.imputer, self.pt, met_data_in, req_sites_dewtemp, reference_sites,
                                                 stations, date_index,reference_station_number, station_list,'dewtemp')

        return(met_data_out_temp,met_data_out_pressure,met_data_out_dewtemp)


    def sort_datasets(self, orig_data_in,req_sites_list,date_index,var_string):

        # add the Date index
        indexed_orig_data = orig_data_in.set_index('Date')

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
            ts = ts.reset_index().set_index(['Date','siteID'])

            # copy data to new array
            full_data_out = full_data_out.append(ts)

        return(full_data_out)


    def organise_data(self, met_data_in,start_date_in,end_date_in,\
                                     req_sites_temp, req_sites_pres, req_sites_dewtemp):

        date_index = pd.date_range(start=start_date, end=end_date, freq='1H', name='Date')

        print('sorting temperature data')
        met_data_out_temp = sort_datasets(met_data_in,req_sites_temp,date_index,'temperature')

        print('sorting pressure data')
        met_data_out_pressure = sort_datasets(met_data_in,req_sites_pres,date_index,'pressure')


        print('sorting dew point temperature')
        met_data_out_dewtemp = sort_datasets(met_data_in,req_sites_dewtemp,date_index,'dewtemp')



        return(met_data_out_temp,met_data_out_pressure,met_data_out_dewtemp)



    def remove_low_temperature_data(self, met_data_in,min_temperature):

        md_cold = met_data_in[met_data_in['temperature']<min_temperature]

        if len(md_cold) > 0:
            print('     deleting this temperature and dew point temperature data')
            print(md_cold)
            met_data_in.loc[md_cold.index,['temperature','dewtemp']] = np.nan
        else:
            print('    no low temperature data to remove')

        return(met_data_in)




    #%%

    def extract_mean_max(self, ts_in,var_in_string,var_out_string):

        tempgroups = ts_in.groupby([pd.Grouper(level="Date",freq='1D'),'siteID'])

        out_data = pd.DataFrame()

        out_data['{}.max'.format(var_out_string)] = tempgroups.max()[var_in_string]
        out_data['{}.mean'.format(var_out_string)] = tempgroups.mean()[var_in_string]
        out_data['{}.flag'.format(var_out_string)] = tempgroups.mean()['{}.flag'.format(var_in_string)]

        return(out_data)


    def combine_and_organise_mean_max(self, met_data_temp,met_data_pres,met_data_rh):

        met_groups_rh   = extract_mean_max(met_data_rh,'rh','RelativeHumidity')
        met_groups_temp = extract_mean_max(met_data_temp,'temperature','Temperature')
        met_groups_pres = extract_mean_max(met_data_pres,'pressure','Pressure')

        combined_data = met_groups_temp.merge(met_groups_rh,how='outer',left_index=True,right_index=True)
        combined_data = combined_data.merge(met_groups_pres,how='outer',left_index=True,right_index=True)

        combined_data.sort_index(level=1,inplace=True)
        combined_data.index = combined_data.index.set_levels(['{} [WEATHER]'.format(x) for x in combined_data.index.levels[1]], level=1)
        combined_data.index.rename(['time_stamp','sensor_name'],inplace=True)

        return(combined_data)


#%%


if __name__ == '__main__':
    # Todo add arguments/processing

    impute_data = True
    print_stats = True
    file_in = 'data_met/temp_rh_press_dewtemp_2016-2019.csv'
    out_dir = 'met_postprocessing'
    file_out = 'daily_mean_max_temp_RH_pres.csv'
    station_drop_list = [117]
    min_temperature = -20
    date_range = [np.datetime64('2016-01-01 00'), np.datetime64('2019-12-31 23')]
    reference_station_number = 5
    verbose = 0

    post_processor = MetPostProcessor(out_dir, verbose)
    post_processor.post_process(file_in, date_range, file_out,  station_drop_list, min_temperature,
                                reference_station_number, impute_data, print_stats)