#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:59:59 2020

Script for testing the imputation of data created using the met data processing script.

This will load data and routines from that script (need to move them )

@author: mbessdl2
"""

import pandas as pd
import numpy as np
from datetime import datetime

import metpy.calc as mpcalc
from metpy.units import units

import geopy.distance as distance

from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge
from sklearn import preprocessing


from sklearn.metrics import mean_squared_error as mse
import seaborn as sns
sns.set_theme()

#%%

def calculate_slices_to_lose(site_indexes,data_loss_fraction):
    
    first_index = len(site_indexes)*data_loss_fraction

    lost_indexes = site_indexes[0:int(first_index)]
    
    return(lost_indexes)


def delete_data_periods(met_data_in,req_sites,data_loss_fraction,var_trim_list):

    met_data_out = met_data_in.copy()
    
    for site in req_sites:
        site_indexes = met_data_in[met_data_in['siteID']==site].index
        lost_indexes = calculate_slices_to_lose(site_indexes,data_loss_fraction)
        
        met_data_out.loc[lost_indexes,var_trim_list] = np.nan
        
    
    return(met_data_out)




#%% obtain stats for each site
    
def stat_check(met_data,met_test,req_sites,var_string,flag_string):
    
    print('stats for species: {}'.format(var_string))
    
    for site in req_sites:
        print('for site {}'.format(site))
        # select site data
        met_data_a = met_data.loc[(slice(None),site),var_string]
        met_test_a = met_test.loc[(slice(None),site),var_string]
        
        met_flag_a = met_test.loc[(slice(None),site),flag_string]
        
        # keep only data which has been imputed
        met_data_a = met_data_a[met_flag_a==1]
        met_test_a = met_test_a[met_flag_a==1]
        
        # remove datapoints that where NaN in original data
        met_test_a = met_test_a[met_data_a.notna()]
        met_data_a = met_data_a[met_data_a.notna()]
        
        
        # calculate the Mean Square Error
        mserror = mse(met_data_a,met_test_a)
        print('  MSError is {}'.format(mserror))
        
        # plot scatter
        met_combined = pd.DataFrame()
        met_combined[var_string] = met_data_a
        met_combined['{}.imputed'.format(var_string)] = met_test_a
        
        sns.jointplot(data=met_combined, x=var_string, y='{}.imputed'.format(var_string))
        
        
    


#%%



if __name__ == '__main__':
    
    file_in = 'data_met/temp_rh_press_dewtemp_2016-2019.csv'
    file_out = 'daily_mean_max_temp_RH_pres.csv'
    col_list = ['date','siteID','temperature','rh','pressure','dewtemp']
    station_drop_list = [117]
    min_temperature = -20
    print_stats = True
    impute_data = True

    data_loss_fraction = 0.5
    var_trim_list = ['temperature','pressure','dewtemp']

    stations = pd.read_csv("station_data/station_data_clean.csv")
    stations = stations.set_index('Station')


    start_date = np.datetime64('2016-01-01 00')
    end_date   = np.datetime64('2019-12-31 23')

    reference_station_number = 5

    imputer = IterativeImputer(random_state=0, add_indicator=True, 
                            initial_strategy='mean',max_iter=300, verbose=0,
                            estimator=BayesianRidge())

    # set the power transform options
    pt = preprocessing.QuantileTransformer(output_distribution='normal',random_state=0)

    
    
    print('loading met data file')
    met_data = load_met_data(file_in,col_list)
    
    print('dropping duplicate values and unwanted stations')
    met_data = find_and_drop_duplicates_and_unwanted_stations(met_data,station_drop_list)
    
    print('dropping single daily measurement stations')
    met_data = drop_single_daily_measurement_stations(met_data,print_stats)
    
    print('filtering to remove unrealistically low temperatures')
    met_data = remove_low_temperature_data(met_data,min_temperature)

    print('filter for minimum data lengths, and reduce dataset to only stations of interest')
    met_data, reference_sites, req_sites_temp, req_sites_pres, req_sites_dewtemp = list_required_and_reference_sites(met_data)


    print('selecting stations to test')
    # replace this with a random number generator?
    req_sites_all = reference_sites[3:40:5]
    reduced_ref_sites = [x for x in reference_sites if x not in req_sites_all]
    

    print('removing defined sections of data for these stations')
    met_data_test = delete_data_periods(met_data,req_sites_all,data_loss_fraction,var_trim_list)


    print('sorting original data (no imputation), returning hourly data')
    met_data_temp,met_data_pres,met_data_dewtemp = organise_data(met_data,start_date,end_date,\
                                    req_sites_all, req_sites_all, req_sites_all)    
        
    print('imputation of the reduced data set, returning hourly data')
    met_test_temp,met_test_pres,met_test_dewtemp = organise_data_imputation(met_data_test,start_date,end_date,stations,\
                                    reduced_ref_sites, req_sites_all, req_sites_all, req_sites_all, \
                                    reference_station_number,imputer,pt)    
    


    
    print('calculation of relative humidity from temperature and dew point temperature')
    met_data_rh = rh_calculations(met_data_temp,met_data_dewtemp,print_stats,met_data)
    met_test_rh = rh_calculations(met_test_temp,met_test_dewtemp,print_stats,met_data)
    
    # calculate the daily max and mean for each station
    met_data_hourly = combine_and_organise_mean_max(met_data_temp,met_data_pres,met_data_rh)
    met_test_hourly = combine_and_organise_mean_max(met_test_temp,met_test_pres,met_test_rh)


    for station in req_sites_all[0:1]:
        stat_check(met_data_temp,met_test_temp,[station],'temperature','temperature.flag')
        stat_check(met_data_pres,met_test_pres,[station],'pressure','pressure.flag')
        stat_check(met_data_dewtemp,met_test_dewtemp,[station],'dewtemp','dewtemp.flag')
        stat_check(met_data_rh,met_test_rh,[station],'rh','rh.flag')

        stat_check(met_data_hourly,met_test_hourly,['{} [WEATHER]'.format(station)],'Temperature.mean','Temperature.flag')
        stat_check(met_data_hourly,met_test_hourly,['{} [WEATHER]'.format(station)],'Temperature.max','Temperature.flag')
        stat_check(met_data_hourly,met_test_hourly,['{} [WEATHER]'.format(station)],'Pressure.mean','Pressure.flag')
        stat_check(met_data_hourly,met_test_hourly,['{} [WEATHER]'.format(station)],'Pressure.max','Pressure.flag')
        stat_check(met_data_hourly,met_test_hourly,['{} [WEATHER]'.format(station)],'RelativeHumidity.mean','RelativeHumidity.flag')
        stat_check(met_data_hourly,met_test_hourly,['{} [WEATHER]'.format(station)],'RelativeHumidity.max','RelativeHumidity.flag')








