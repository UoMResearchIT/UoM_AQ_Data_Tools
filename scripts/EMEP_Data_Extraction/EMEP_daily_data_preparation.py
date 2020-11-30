#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:11:54 2020

@author: mbessdl2
"""

from pathlib import Path
import pandas as pd



#%%

def load_tidy_emep_data(emep_file,spc_list):
    
    emep_dataframe = pd.read_csv(emep_file)

    # correct name of NOx variable
    emep_dataframe = emep_dataframe.rename(columns={'NOx':'NOXasNO2'})
    
    # convert Date string to datetime object, and set as index
    emep_dataframe['Date'] = pd.to_datetime(emep_dataframe['Date'])
    emep_dataframe = emep_dataframe.set_index('Date')

    # change the column headers for the chemical species
    emep_dataframe = emep_dataframe.rename(columns=spc_dict)

    return(emep_dataframe)


def gather_and_organise_emep_data(emep_dataframe_internal,site_list_internal,spc_list):
 
    final_dataframe = pd.DataFrame()

    for site in site_list_internal:
        # select our subset of metadata for this station
        print("processing site {}".format(site))

        working_hourly_dataframe = emep_dataframe_internal[emep_dataframe_internal["SiteID"] == site]

        # postprocessing the data set, to get daily data
        daily_grouped_data = working_hourly_dataframe.groupby(pd.Grouper(freq='1D'))
        temp_dataframe = pd.DataFrame()
        for spc in spc_list:
            temp_dataframe['{}.mean'.format(spc)] = daily_grouped_data.mean()[spc]
            temp_dataframe['{}.max'.format(spc)] = daily_grouped_data.max()[spc]
            temp_dataframe['sensor_name'] = '{} [AQ]'.format(site)
        
        temp_dataframe = temp_dataframe.rename_axis('time_stamp',axis=0).set_index('sensor_name',append=True)
        final_dataframe = final_dataframe.append(temp_dataframe)


    return(final_dataframe)



if __name__ == '__main__':



    # set the species list
    spc_list = ['O3','NO2','SO2','NOXasNO2','PM2.5','PM10']
    # create the EMEP variable names from these
    spc_emep_list = list('{}_EMEP'.format(spc) for spc in spc_list)
    spc_dict = dict((spc,'{}_EMEP'.format(spc)) for spc in spc_list)


    # emep data file
    emep_file = Path("emep_site_data_hourly_2016-2019.csv")

    # load emep data, rename the NOX columns, and convert Date strings to datetime objects
    emep_dataframe = load_tidy_emep_data(emep_file,spc_dict)

    # obtain the list of sites
    site_list = emep_dataframe['SiteID'].unique()

    # create our daily max and mean dataset
    daily_dataframe = gather_and_organise_emep_data(emep_dataframe,site_list,spc_emep_list)

    # write to file
    daily_dataframe.to_csv(Path('emep_daily_data_2016-2019.csv'),index=True,header=True,float_format='%.2f')



