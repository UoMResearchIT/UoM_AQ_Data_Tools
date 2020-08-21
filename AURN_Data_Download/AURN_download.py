#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script for automated downloading of AURN data for a given date range.

"""

import wget
from pathlib import Path
import os
import os.path
import pyreadr
import datetime
import pandas as pd
import numpy as np



# functions for downloading data

def define_years_to_download(subset_df,years):

    # Check to see if your requested years will work and if not, change it
    # to do this lets create two new columns of datetimes for earliest and latest
    datetime_start=pd.to_datetime(subset_df['start_date'].values, format='%Y/%m/%d').year
    #Problem with the end date is it could be ongoing. In which case, convert that entry into a date and to_datetime
    now = datetime.datetime.now()
    datetime_end_temp=subset_df['end_date'].values
    step=0
    for i in datetime_end_temp:
        if i == 'ongoing':
            datetime_end_temp[step]=str(now.year)+'-'+str(now.month)+'-'+str(now.day)
        step+=1
    datetime_end = pd.to_datetime(datetime_end_temp).year

    earliest_year = np.min(datetime_start)
    latest_year = np.max(datetime_end)
    
    # now create list of years to process
    years_process = []
    
    for year in years:
        if(year <= latest_year and year >= earliest_year):
            years_process.append(year)
    
    return(years_process)


def download_and_open_datafiles(subset_df,site,station_name,years,data_path):
    
    downloaded_site_data = []

    for year in years:
        try:
            downloaded_file=site+"_"+str(year)+".RData"
            download_url = "https://uk-air.defra.gov.uk/openair/R_data/"+downloaded_file
            print("\tdownloading file {}".format(download_url))

            # Check to see if file exists or not. Special case for current year as updates on hourly basis
            filename_path = data_path.joinpath(downloaded_file)
            if (filename_path.is_file() is True):
                print("\t\tData file already exists, will use this")
            else:
                print("\t\tDownloading data file for ", station_name ," in ",str(year))
                wget.download(download_url,out=str(data_path))

            # Read the RData file into a Pandas dataframe
            downloaded_data = pyreadr.read_r(str(filename_path))
            # Drop non-required fields
            


            # Append to dataframe list
            downloaded_site_data.append(downloaded_data[site+"_"+str(year)])
        except:
            print("\t\tCouldn't download and extract data from {} for {}".format(year,station_name))


    return(downloaded_site_data)


def load_sort_data(downloaded_site_data):

    final_dataframe = pd.concat(downloaded_site_data, axis=0, ignore_index=True)

    final_dataframe['datetime'] = pd.to_datetime(final_dataframe['date'])
    final_dataframe=final_dataframe.sort_values(by='datetime',ascending=True)
    #final_dataframe=final_dataframe.set_index('datetime')

    return(final_dataframe)



def postprocess_data(hourly_dataframe,site):

    columns_of_interest = ['datetime','O3','NO2','SO2','NOXasNO2','PM2.5','PM10']
    ds_columns = hourly_dataframe.columns
    
    # retain the data we are interested in (as not all datasets have all variables)
    columns_to_retain = set(columns_of_interest) & set(ds_columns)
    working_dataframe = hourly_dataframe[columns_to_retain]
    
    working_dataframe = working_dataframe.rename(columns={'datetime':'Date'})
    tempgroups = working_dataframe.groupby(pd.Grouper(key='Date', freq='1D'))
    
    data_counts = tempgroups.count()
    data_max    = tempgroups.max()
    data_mean   = tempgroups.mean()
    
    cols_old = data_counts.columns
    
    cols_counts = dict((key,key+'_count') for key in cols_old.values) 
    cols_max    = dict((key,key+'_max') for key in cols_old.values) 
    cols_mean   = dict((key,key+'_mean') for key in cols_old.values) 
    
    data_counts = data_counts.rename(columns=cols_counts)
    data_max    = data_max.rename(columns=cols_max)
    data_mean   = data_mean.rename(columns=cols_mean)
    
    data_out = data_mean.join([data_max,data_counts])
    
    # add the site as a new column, and set as part of multiindex with the date
    site_name = "{} [AQ]".format(site)
    
    data_out['siteID'] = site_name
    data_out = data_out.reset_index(drop=False).set_index(['Date','siteID'])    
    
    return(data_out)


def extract_site_data(site_list,metadata,years,data_path,save_to_csv):

    final_dataframe = pd.DataFrame()

    for site in site_list:

        # select our subset of metadata for this station
        subset_df = metadata['AURN_metadata'][metadata['AURN_metadata'].site_id == site]
        station_name = subset_df['site_name'].values[0]

        print("processing site {} ({})".format(site,station_name))
    

        # get the list of years of data for this station
        years_process = define_years_to_download(subset_df,years)

        # if the list of years is empty, then skip this site
        if not years_process:
            print("\t\tNo data for years of interest, skipping this site")
            continue

        # download the datasets
        downloaded_site_data = download_and_open_datafiles(subset_df,site,station_name,years_process,data_path)

        # if we couldn't download the data, skip this site
        if len(downloaded_site_data) == 0:
            print("\t\tNo data could be downloaded for {}".format(station_name))
            continue

        # combine and sort data
        full_hourly_dataframe = load_sort_data(downloaded_site_data)


        # postprocessing the data set, to get daily data
        final_dataframe = final_dataframe.append(postprocess_data(full_hourly_dataframe,site))


        # Now save the full hourly dataframe as a .csv file
        if save_to_csv is True:
            full_hourly_dataframe.to_csv(data_path.joinpath(site+'.csv'), index = False, header=True)



    return(final_dataframe)






if __name__ == '__main__':

    meta_data_url = "https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData"
    meta_data_filename = 'AURN_metadata.RData'

    data_url = "https://uk-air.defra.gov.uk/openair/R_data/"

    # settings and containers for downloading data
    years = [2016,2017,2018,2019]

    # process control flags
    save_to_csv = True


    # Does the metadatafile exist?
    if os.path.isfile(meta_data_filename) is True:
        print("Meta data file already exists in this directory, will use this")
    else:
        print("Downloading Meta data file")
        wget.download(meta_data_url)
    


    # Read the RData file into a Pandas dataframe
    metadata = pyreadr.read_r(meta_data_filename)


    # If a single year is passed then convert to a list with a single value
    if type(years) is int:
        years = [years]
    current_year = datetime.datetime.now()


    base_path = Path("AURN_data_download")
    if (base_path.is_dir() is False): 
        base_path.mkdir() 
    # change this if we need more organised data storage later  
    data_path = base_path

    # get list of sites to process
    site_list = metadata['AURN_metadata']['site_id'].unique()
    
    # pull out the daily mean and max values for the site list
    daily_dataframe = extract_site_data(site_list,metadata,years,data_path,save_to_csv)

    # sort the data
    daily_dataframe = daily_dataframe.sort_index()
    
    # write this dataset to file
    daily_dataframe.to_csv(data_path.joinpath('pollution_daily_data_{}-{}.csv'.format(years[0],years[-1])),index=True,header=True,float_format='%.2f')
















