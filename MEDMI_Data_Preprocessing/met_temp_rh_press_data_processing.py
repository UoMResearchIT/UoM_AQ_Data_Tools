#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 2020

Script for processing the temperature, RH and pressure hourly meteorological
data, to create daily mean and maximum values. Data is selected primarily on 
the availability of temperature data, with relative humidity and pressure being
secondary variables, extracted only where there is also temperature data.

For relative humidity data this does not cause any loss of data compared to
extracting rel hum data as the primary variable. For pressure data there is a
loss in data, but this is small (e.g. only 0.7% of days with 24 hourly data
points, 216,238 rather than 217,784 days).

Dataset is mostly 24 hourly data (see below). Some days and stations are missing 
data (between 1-23 hourly data). There are also some duplicate readings per hour.
Some single readings per day - at 8am (always?) - ~1/3 of the data. These are
synoptic spot readings, and should be removed (as they are no use for calculating
daily statistics). The (vast majority) of the duplicated data are because some
stations take both SYNOP and METAR (aviation) measurements. SYNOP measurements
are more accurate, so we will take these where they are available - and only
use METAR data where (if at all) SYNOP data is not available.

(pers. comms., Martyn Sunter, Met Office, July 2020)

Raw data count per day is summarised below. There are 10,217 day / station
combinations without relative humidity measurements, and 417,559 without 
pressure data. 

total temperature daily data count is: 637195
# data points per day, total daily data point counts
       temp,  rel hum, pressure
   0,      0,  10217, 417559
   1, 179196, 171356,     15
   2,     76,     88,      8
   3,     89,    109,     10
   4,    110,    133,      5
   5,    138,    153,     10
   6,    226,    230,     13
   7,    448,    474,     11
   8,   1082,   1110,     32
   9,   1636,   1636,     37
  10,   1878,   1886,     57
  11,   2657,   2632,     54
  12,   1714,   1612,     47
  13,   1391,   1372,     53
  14,   1691,   1716,     59
  15,   1889,   1921,     67
  16,   1945,   1987,     66
  17,   2132,   2153,     80
  18,    802,    843,    115
  19,    550,    599,    180
  20,    829,    896,    228
  21,   1396,   1464,    360
  22,   2331,   2408,    502
  23,   7425,   7878,   1389
  24, 400900, 398549, 216238
  25,    136,    117,      0
  26,     68,     20,      0
  27,    127,     20,      0
  28,    126,     30,      0
  29,    160,     38,      0
  30,    185,     15,      0
  31,    212,     16,      0
  32,     92,     47,      0
  33,    173,    153,      0
  34,    555,    537,      0
  35,    815,    811,      0
  36,    357,    360,      0
  37,    500,    484,      0
  38,    362,    362,      0
  39,   1070,   1051,      0
  40,    293,    297,      0
  41,    152,    153,      0
  42,    164,    163,      0
  43,    114,    116,      0
  44,     49,     56,      0
  45,     29,     30,      0
  46,     89,    101,      0
  47,    371,    423,      0
  48,  18465,  18373,      0
  

Duplicated hourly readings are filtered out. There are 1,042,218 hourly measurements
which are duplicated. 529,674 of these have no pressure data, while 512,544 do
have pressure data. There are 2 duplicates which have pressure data for both
readings, of these I've retained the 1st value each time. Similarly there are
8,567 duplicates for which neither reading has pressure data. Of these I have
retained the first value each time (as the aim here is (1) to eliminate duplicates,
and (2) retain the higher resolution data if possible; so even where (2) isn't 
possible we make sure that (1) is followed). 

                                    
total temperature daily data count is: 637195
# data points per day, total daily data point counts
       temp, rel hum, pressure
   0,      0,  10218, 417559
   1, 179196, 171356,     15
   2,     76,     88,      8
   3,     89,    109,     10
   4,    110,    133,      5
   5,    138,    153,     10
   6,    226,    231,     13
   7,    449,    475,     11
   8,   1083,   1112,     32
   9,   1636,   1637,     37
  10,   1878,   1887,     57
  11,   2658,   2634,     54
  12,   1715,   1613,     47
  13,   1391,   1374,     53
  14,   1692,   1717,     59
  15,   1890,   1923,     67
  16,   1946,   1987,     66
  17,   2135,   2155,     80
  18,    805,    848,    115
  19,    554,    605,    180
  20,    831,    900,    228
  21,   1402,   1472,    360
  22,   2338,   2429,    503
  23,   7459,   7933,   1389
  24, 425498, 422206, 216237  


Removing days/stations for which temperature readings are 1 reduces the daily
data count to 457,999.  
total temperature daily data count is: 457999
# data points per day, total daily data point counts
       temp, rel hum, pressure
   0,      0,   2359, 238375
   1,      0,     19,      3
   2,     76,     88,      8
   3,     89,    109,     10
   4,    110,    133,      5
   5,    138,    153,     10
   6,    226,    231,     13
   7,    449,    475,     11
   8,   1083,   1112,     32
   9,   1636,   1637,     37
  10,   1878,   1887,     57
  11,   2658,   2634,     54
  12,   1715,   1613,     47
  13,   1391,   1374,     53
  14,   1692,   1717,     59
  15,   1890,   1923,     67
  16,   1946,   1987,     66
  17,   2135,   2155,     80
  18,    805,    848,    115
  19,    554,    605,    180
  20,    831,    900,    228
  21,   1402,   1472,    360
  22,   2338,   2429,    503
  23,   7459,   7933,   1389
  24, 425498, 422206, 216237





@author: mbessdl2
"""
#import appnope
#appnope.nope()


import pandas as pd
from datetime import datetime
import numpy as np

#%% function for creating date objects

def calc_date(date_in):
    return datetime.strptime(date_in,'%Y-%m-%d %H:%M:%S')

def calc_nanmean(data_in):
    return np.nanmean(data_in)

#%% function for writing out some information about the data count stats

def print_data_count_stats(dc_in):

    print('total temperature daily data count is: {}'.format(dc_in.count().values[0]))
    print('# data points per day, total daily data point counts')
    print('                      temperature, rel hum, pressure')
    for dpoint in range(0,49):
        dcount = dc_in[dc_in.temperature==dpoint].count().values[0]
        ecount = dc_in[dc_in.rh==dpoint].count().values[0]
        fcount = dc_in[dc_in.pressure==dpoint].count().values[0]
        print('{0:4},{1:7},{2:7},{3:7}'.format(dpoint,dcount,ecount,fcount))




#%% read in data

#file_in = 'data_met/temp_rh_press_2016-2019.csv'
file_in = 'data_met/temp.data.csv'


met_data = pd.read_csv(file_in,usecols=['date','siteID','temperature','rh','pressure'],skiprows=1)


#%% create proper date objects, and drop old data strings

met_data['Date'] = met_data['date'].apply(calc_date)

met_data = met_data.drop(columns='date')


#%% find duplicate indexes

#met_data = met_data.set_index(['Date','siteID'])

met_duplicates = met_data[met_data.duplicated(subset=['Date','siteID'],keep=False)]

met_dup_with_pres = met_duplicates[met_duplicates['pressure'].isna()==False]
met_dup_no_pres = met_duplicates[met_duplicates['pressure'].isna()==True]

met_dup_with_pres_dups = met_dup_with_pres[met_dup_with_pres.duplicated(subset=['Date','siteID'],keep='first')]

met_dup_no_pres_dups = met_dup_no_pres[met_dup_no_pres.duplicated(subset=['Date','siteID'],keep='last')]
met_dup_no_pres_reduced = met_dup_no_pres.drop(index=met_dup_no_pres_dups.index)


#%% drop the unwanted indexes from our original met_data
#   this eliminates all duplicated data points

indexes_to_drop = met_dup_with_pres_dups.index
indexes_to_drop = indexes_to_drop.append(met_dup_no_pres_reduced.index)

met_data_reduced = met_data.drop(index=indexes_to_drop)



#%% group by station ID and date (at a 1 day frequency)
#   and get some data counts for these groups

tempgroups = met_data_reduced.groupby(['siteID',pd.Grouper(key='Date', freq='1D')])

data_counts = tempgroups.count()

#%% print out some stats about the data point counts for each day

print_data_count_stats(data_counts)

    
#%% drop the datasets which we're not interested in.
#   1) drop all datasets with a single daily value

trimmed_data_counts = data_counts[data_counts['temperature']>1]
low_data_counts = trimmed_data_counts[trimmed_data_counts['temperature']<24]

full_day_counts = data_counts[data_counts['temperature']==24]

halfplus_day_counts = data_counts[data_counts['temperature']>=12]
threequarterplus_day_counts = data_counts[data_counts['temperature']>=18]
twothirdsplus_day_counts = data_counts[data_counts['temperature']>=16]


trimmed_data_counts_rh = data_counts[data_counts['rh']>1]
full_day_counts_rh = data_counts[data_counts['rh']==24]
halfplus_day_counts_rh = data_counts[data_counts['rh']>=12]

trimmed_data_counts_pres = data_counts[data_counts['pressure']>1]
full_day_counts_pres = data_counts[data_counts['pressure']==24]
halfplus_day_counts_pres = data_counts[data_counts['pressure']>=12]


###### switch over to MICE_temp_rh_pres.py script

#%% create database for recording how many of each data count for each site
#  index: site
#  columns: counts per day (1 to 48)

site_data_count = pd.DataFrame(index=data_counts.index.get_level_values(0).unique(),columns=list(range(1,49)))

site_data_count.iloc[:,:] = 0

#%% calculate the counts for each site

for site in site_data_count.index:
    data_temp = data_counts.loc[(site,)]
    for dpoint in range(1,49):
        dcount = data_temp[data_temp.temperature==dpoint].count()
        site_data_count.loc[site,dpoint] = dcount.values[0]

site_data_count.loc[:,'none 1']=site_data_count.iloc[:,1:48].sum(axis=1)


site_data_count[site_data_count['none 1']==0].to_csv('sites_with_single_daily_values.csv')
site_data_count[site_data_count['none 1']>0].to_csv('sites_with_multiple_daily_values.csv')




#%% calculate the daily mean and max values, and data count

temp_daily_data = temperature_data_hourly.groupby(['siteID',pd.Grouper(key='Date', freq='1D0H',base=8,loffset='-8H')]).agg(['mean','max','count'])
temp_daily_data = temperature_data.groupby(['siteID',pd.Grouper(key='Date', freq='24H',base=8,loffset='-8H')]).agg(['mean','max','count'])




