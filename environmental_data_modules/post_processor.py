#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from datetime import datetime
from abc import ABCMeta, abstractmethod
import os

try:
    import pandas as pd
    import numpy as np
    import geopy.distance as distance
except:
    pass

from environmental_data_modules import EnvironmentModule

class PostProcessor(EnvironmentModule):
    __metaclass__ = ABCMeta

    DEFAULT_IMPUTE_DATA = True
    DEFAULT_PRINT_STATS = True
    DEFAULT_SKIP_INPUT_ROWS = 1
    DEFAULT_MIN_YEARS = 1


    def __init__(self, out_dir=EnvironmentModule.DEFAULT_OUT_DIR, verbose=EnvironmentModule.DEFAULT_VERBOSE):
        super(PostProcessor, self).__init__(out_dir, verbose)

        self._impute_data = PostProcessor.DEFAULT_IMPUTE_DATA
        self._print_stats = PostProcessor.DEFAULT_PRINT_STATS
        self._file_in = None
        self._skip_input_rows = PostProcessor.DEFAULT_SKIP_INPUT_ROWS
        self._min_years = PostProcessor.DEFAULT_MIN_YEARS
        self._stations = None


    @property
    def stations(self):
        return self._stations

    @property
    def imputer(self):
        return self._imputer

    @imputer.setter
    def imputer(self, imputer):
        if imputer is None or type(imputer).__name__ == 'IterativeImputer':
            self._imputer = imputer
        else:
            raise ValueError('Error setting imputer, incorrect object type: {}'.format(type(imputer).__name__))

    @property
    def transformer(self):
        return self._transformer

    #%% station geographic routines
    def calc_station_distances(self, stations_in, stat_location):

        station_distances = pd.DataFrame(index=stations_in.index)
        station_distances['Distance'] = np.nan

        for index, row in stations_in.iterrows():
            new_location = (row['Latitude'],row['Longitude'])

            station_distances.loc[index]['Distance'] = distance.distance(stat_location,new_location).km

        return station_distances


    # %% function for creating date objects
    def parse_calcs_date(self, date_in):
        return datetime.strptime(date_in, self._date_calcs_format)

    def calc_nanmean(self, data_in):
        return np.nanmean(data_in)
