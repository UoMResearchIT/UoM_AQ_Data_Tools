#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:40:42 2020

wind data - it's looking reasonable?!?!

total wind data count is: 233910
# data points per day, total data point count
   1,     28
   2,     22
   3,     34
   4,     24
   5,     31
   6,     30
   7,     50
   8,     70
   9,     81
  10,    112
  11,    112
  12,    119
  13,     87
  14,    110
  15,    124
  16,    117
  17,     86
  18,    121
  19,    157
  20,    187
  21,    273
  22,    442
  23,   1146
  24, 230347

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

file_in = 'data_met/wind_2016-2019.csv'

wind_data = pd.read_csv(file_in,usecols=['date','siteID','windspeed','winddir'],skiprows=1)

#%% create proper date objects, and drop old data strings

wind_data['Date'] = wind_data['date'].apply(calc_date)

wind_data = wind_data.drop(columns='date')

#%% group by station ID and date (at a 1 day frequency)
#   and get some data counts for these groups

windgroups = wind_data.groupby(['siteID',pd.Grouper(key='Date', freq='1D')])

data_counts = windgroups.count()

#%% print out some stats about the data point counts for each day

print('total wind data count is: {}'.format(data_counts.count().values[0]))
print('# data points per day, total data point count')
for dpoint in range(1,25):
    dcount = data_counts[data_counts.winddir==dpoint].count()
    print('{0:4},{1:7}'.format(dpoint,dcount.values[0]))

