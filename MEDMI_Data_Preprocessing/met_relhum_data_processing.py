#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 16:40:42 2020

total relative humidity data count is: 626978
# data points per day, total data point count
   1, 171356
   2,     88
   3,    109
   4,    133
   5,    153
   6,    230
   7,    474
   8,   1110
   9,   1636
  10,   1886
  11,   2632
  12,   1612
  13,   1372
  14,   1716
  15,   1921
  16,   1987
  17,   2153
  18,    843
  19,    599
  20,    896
  21,   1464
  22,   2408
  23,   7878
  24, 398549
  25,    117
  26,     20
  27,     20
  28,     30
  29,     38
  30,     15
  31,     16
  32,     47
  33,    153
  34,    537
  35,    811
  36,    360
  37,    484
  38,    362
  39,   1051
  40,    297
  41,    153
  42,    163
  43,    116
  44,     56
  45,     30
  46,    101
  47,    423
  48,  18373

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

file_in = 'data_met/relhum_2016-2019.csv'

relhum_data = pd.read_csv(file_in,usecols=['date','siteID','rh'],skiprows=1)

#%% create proper date objects, and drop old data strings

relhum_data['Date'] = relhum_data['date'].apply(calc_date)

relhum_data = relhum_data.drop(columns='date')

#%% group by station ID and date (at a 1 day frequency)
#   and get some data counts for these groups

rhgroups = relhum_data.groupby(['siteID',pd.Grouper(key='Date', freq='1D')])

data_counts = rhgroups.count()

#%% print out some stats about the data point counts for each day

print('total relative humidity data count is: {}'.format(data_counts.count().values[0]))
print('# data points per day, total data point count')
for dpoint in range(1,100):
    dcount = data_counts[data_counts.rh==dpoint].count()
    print('{0:4},{1:7}'.format(dpoint,dcount.values[0]))

