#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:45:13 2020

Script for processing the temperature hourly meteorological data, 
to create daily mean and maximum values.


Dataset is mostly 24 hourly data (see below). Some missing data (between 1-23 hourly data).

Some duplicate readings per hour - unclear what these are.

Some single readings per day - at 8am (always?) - ~1/3 of the data.


total temperature daily data count is: 637195
# data points per day, total daily data point count
   1, 179196
   2,     76
   3,     89
   4,    110
   5,    138
   6,    226
   7,    448
   8,   1082
   9,   1636
  10,   1878
  11,   2657
  12,   1714
  13,   1391
  14,   1691
  15,   1889
  16,   1945
  17,   2132
  18,    802
  19,    550
  20,    829
  21,   1396
  22,   2331
  23,   7425
  24, 400900
  25,    136
  26,     68
  27,    127
  28,    126
  29,    160
  30,    185
  31,    212
  32,     92
  33,    173
  34,    555
  35,    815
  36,    357
  37,    500
  38,    362
  39,   1070
  40,    293
  41,    152
  42,    164
  43,    114
  44,     49
  45,     29
  46,     89
  47,    371
  48,  18465

@author: mbessdl2
"""


import pandas as pd
from datetime import datetime
import numpy as np

#%% function for creating date objects

def calc_date(date_in):
    return datetime.strptime(date_in,'%Y-%m-%d %H:%M:%S')

def calc_nanmean(data_in):
    return np.nanmean(data_in)


#%% read in data

file_in = 'data_met/temperature_2016-2019.csv'

temperature_data = pd.read_csv(file_in,usecols=['date','siteID','temperature'],skiprows=1)


#%% create proper date objects, and drop old data strings

temperature_data['Date'] = temperature_data['date'].apply(calc_date)

temperature_data = temperature_data.drop(columns='date')

#%% group by station ID and date (at a 1 day frequency)
#   and get some data counts for these groups

tempgroups = temperature_data.groupby(['siteID',pd.Grouper(key='Date', freq='1D')])

data_counts = tempgroups.count()

#%% print out some stats about the data point counts for each day

print('total temperature daily data count is: {}'.format(data_counts.count().values[0]))
print('# data points per day, total daily data point count')
for dpoint in range(1,49):
    dcount = data_counts[data_counts.temperature==dpoint].count()
    print('{0:4},{1:7}'.format(dpoint,dcount.values[0]))

#%% 
