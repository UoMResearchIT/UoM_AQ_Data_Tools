from abc import ABCMeta, abstractmethod
from cmath import polar
import argparse
import json
import os, errno
from datetime import datetime


class EnvironmentWorkflow:
    __metaclass__ = ABCMeta
    DEFAULT_OUT_FILE_SUFFIX = ''
    DEFAULT_OUT_DIR = 'met_extracted_data'
    DEFAULT_DATE_RANGE = ['2016-1-1 0', '2019-12-31 23']
    UK_LATITUDES = [48., 60.]
    UK_LONGITUDES = [-11., 3.]
    DEFAULT_VERBOSE = 0
    DEFAULT_DATE_RANGE_FORMAT = '%Y-%m-%d %H'
    LONGITUDE_RANGE = [-180., 360.]
    LATITUDE_RANGE = [-90., 90.]
    DEFAULT_COLS_BASE = ['date', 'siteID']


    def __init__(self, dir_name=DEFAULT_OUT_DIR, verbose=DEFAULT_VERBOSE):
        self.out_dir = dir_name
        self.verbose = verbose
        self.latitude_range = EnvironmentWorkflow.UK_LATITUDES
        self.longitude_range = EnvironmentWorkflow.UK_LONGITUDES
        self._date_range_format = EnvironmentWorkflow.DEFAULT_DATE_RANGE_FORMAT
        self.date_range = EnvironmentWorkflow.DEFAULT_DATE_RANGE
        self._file_out = None
        self._cols_base = EnvironmentWorkflow.DEFAULT_COLS_BASE
        self._cols_specific = []


    @staticmethod
    def all_subclasses(cls):
        return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in EnvironmentWorkflow.all_subclasses(c)])

    @property
    def date_range(self):
        return self.__date_range

    @date_range.setter
    def date_range(self, range):
        # Example input: ['2017-1-1 0', '2017-06-30 23']
        try:
            datetime_1 = datetime.strptime(range[0], self._date_range_format)
            datetime_2 = datetime.strptime(range[1], self._date_range_format)
        except ValueError as err:
            raise err
        if datetime_1 >= datetime_2:
            raise ValueError('Start date is not earlier than end date.')
        self.__date_range = [range[0], range[1]]

    @property
    def latitude_range(self):
        return self.__latitude_range

    @latitude_range.setter
    def latitude_range(self, range):
        try:
            val_1 = float(range[0])
            val_2 = float(range[1])
        except ValueError as err:
            raise err
        if not ((EnvironmentWorkflow.LATITUDE_RANGE[0] <= val_1) and (val_1 <= EnvironmentWorkflow.LATITUDE_RANGE[1])):
            raise ValueError('Latitude first value falls outside global range')
        if not ((EnvironmentWorkflow.LATITUDE_RANGE[0] <= val_2) and (val_2 <= EnvironmentWorkflow.LATITUDE_RANGE[1])):
            raise ValueError('Latitude last value falls outside global range')
        self.__latitude_range = [val_1, val_2]

    @property
    def longitude_range(self):
        return self.__longitude_range

    @longitude_range.setter
    def longitude_range(self, range):
        try:
            val_1 = float(range[0])
            val_2 = float(range[1])
        except ValueError as err:
            raise err
        if not (EnvironmentWorkflow.LONGITUDE_RANGE[0] <= val_1 <= EnvironmentWorkflow.LONGITUDE_RANGE[1]):
            raise ValueError('Longitude first value falls outside global range')
        if not (EnvironmentWorkflow.LONGITUDE_RANGE[0] <= val_2 <= EnvironmentWorkflow.LONGITUDE_RANGE[1]):
            raise ValueError('Longitude last value falls outside global range')
        self.__longitude_range = [val_1, val_2]

    @property
    def out_dir(self):
        return self.__out_dir

    @out_dir.setter
    def out_dir(self, dir_name):
        try:
            dir_name = str(dir_name)
        except ValueError as err:
            raise err
        try:
            os.makedirs(dir_name)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise ValueError("Directory name {} cannot be created.".format(dir_name))
        self.__out_dir = dir_name

    @property
    def verbose(self):
        return self.__verbose

    @verbose.setter
    def verbose(self, verbose):
        try:
            verbose = int(verbose)
        except ValueError as err:
            raise err
        self.__verbose = max(0, verbose)

    @property
    def file_out(self):
        return self._file_out

    def set_region(self, latitude_range, longitude_range):
        self.latitude_range = latitude_range
        self.longitude_range = longitude_range

    def get_all_cols(self):
        return self._cols_base + self._cols_specific