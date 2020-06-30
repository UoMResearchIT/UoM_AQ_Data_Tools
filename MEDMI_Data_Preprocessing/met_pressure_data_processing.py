#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:40:42 2020

total pressure data count is: 220307
# data points per day, total data point count
   1,     13
   2,      8
   3,      9
   4,      5
   5,      9
   6,     11
   7,      9
   8,     29
   9,     34
  10,     49
  11,     51
  12,     44
  13,     52
  14,     51
  15,     58
  16,     60
  17,     67
  18,     99
  19,    147
  20,    157
  21,    208
  22,    288
  23,   1065
  24, 217784

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

file_in = 'data_met/stnpres_2016-2019.csv'

pressure_data = pd.read_csv(file_in,usecols=['date','siteID','pressure'],skiprows=1)

#%% create proper date objects, and drop old data strings

pressure_data['Date'] = pressure_data['date'].apply(calc_date)

pressure_data = pressure_data.drop(columns='date')

#%% group by station ID and date (at a 1 day frequency)
#   and get some data counts for these groups

presgroups = pressure_data.groupby(['siteID',pd.Grouper(key='Date', freq='1D')])

data_counts = presgroups.count()

#%% print out some stats about the data point counts for each day

print('total pressure data count is: {}'.format(data_counts.count().values[0]))
print('# data points per day, total data point count')
for dpoint in range(1,100):
    dcount = data_counts[data_counts.pressure==dpoint].count()
    print('{0:4},{1:7}'.format(dpoint,dcount.values[0]))

