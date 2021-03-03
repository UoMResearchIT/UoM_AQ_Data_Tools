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
    """
        Abstract class, parent of classes used for post-processing data that has been extracted from AURN/MEDMI servers.
    """
    __metaclass__ = ABCMeta

    DEFAULT_IMPUTE_DATA = True
    DEFAULT_PRINT_STATS = True
    DEFAULT_SKIP_INPUT_ROWS = 1
    DEFAULT_SAVE_TO_CSV = True


    def __init__(self, out_dir=EnvironmentModule.DEFAULT_OUT_DIR, verbose=EnvironmentModule.DEFAULT_VERBOSE):
        """ Initialise instance of the PostProcessor class.
            Initialises the private class variables

            Args:
                out_dir: (string) directory to be used for all outputs
                verbose: (integer) level of verbosity in output.

            Returns:
                Initialised instance of subclass of PostProcessor

        """
        super(PostProcessor, self).__init__(out_dir, verbose)

        if type(self) == PostProcessor:
            raise NotImplementedError('PostProcessor cannot be instantiated directly')

        self._impute_data = PostProcessor.DEFAULT_IMPUTE_DATA
        self._print_stats = PostProcessor.DEFAULT_PRINT_STATS
        self._file_in = None
        self._skip_input_rows = PostProcessor.DEFAULT_SKIP_INPUT_ROWS
        self._station_data = None
        self._min_years_reference = None
        self._min_years = None
        self._imputer = None

    @abstractmethod
    def process(self, file_in, outfile_suffix='', date_range=None, impute_data=DEFAULT_IMPUTE_DATA,
                save_to_csv=DEFAULT_SAVE_TO_CSV):
        raise NotImplementedError("Must override process")

    @abstractmethod
    def impute_method_setup(self):
        raise NotImplementedError("Must override impute_method_setup")

    ### Class properties: Get/Sets ###

    @property
    def station_data(self):
        return self._station_data

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

    @property
    def impute_data(self):
        return self._impute_data

    @impute_data.setter
    def impute_data(self, impute):
        if not impute in [True, False]:
            raise Exception('impute_data value must be a boolean. {} was input.'.format(impute))
        if impute and self.imputer is None:
            raise Warning('Imputation requested but imputer is currently None. \n Please run impute_method_setup')
        self._impute_data = impute

    @property
    def min_years_reference(self):
        return self._min_years_reference

    @min_years_reference.setter
    def min_years_reference(self, min_years):
        try:
            min_years = float(min_years)
        except:
            raise ValueError('min_years_reference value ({}) must be numeric'.format(min_years))
        if min_years < 0:
            raise ValueError('min_years_reference must be non-negative.')
        self._min_years_reference = min_years

    @property
    def min_years(self):
     return self._min_years

    @min_years.setter
    def min_years(self, min_years):
        try:
            min_years = float(min_years)
        except:
            raise ValueError('min_years value ({}) must be numeric'.format(min_years))
        if min_years < 0:
            raise ValueError('min_years must be non-negative.')
        self._min_years = min_years

    ### station geographic routines

    def calc_station_distances(self, stations_in, location):
        """
        Calculates the distances between a list of stations and a specified location.
        
        Args:
            stations_in: pandas.Dataframe containing station locations
                Required Index:
                    site_id   (str): identifiers for the stations
                Required Columns:
                    latitude (float):
                    longitude (float):
            location (tuple, float): (latitude, longitude) of station of interest
        
        Returns:
            station_distances: pandas.Dataframe containing (sorted) list of distances to location for each station
                Required Index:
                    site_id   (str): identifiers for the stations
                Required Columns:
                    distance (float): distance to listed station from station of interest, in km
        """
        if stations_in is None:
            raise ValueError("stations_in must not be None")
        if location is None:
            raise ValueError("location must not be None")
        if not isinstance(stations_in, pd.DataFrame):
            raise ValueError('stations_in must be a valid pandas Dataframe')
        if not isinstance(location, tuple):
            raise ValueError('location must be a valid tuple')
        if len(location) != 2:
            raise ValueError('location must contain 2 items')
        if not all(isinstance(item, float) for item in location):
            raise ValueError('both elements in location list must be floats')
        if stations_in.index.dtype != str:
            ValueError('stations_in dataframe must have index with string type (containing site_id\'s')

        station_distances = pd.DataFrame(index=stations_in.index)
        station_distances.index.names = ['site_id']
        station_distances['distance'] = np.nan

        for index, row in stations_in.iterrows():
            new_location = (row['latitude'], row['longitude'])
            station_distances.loc[index]['distance'] = distance.distance(location, new_location).km

        return station_distances


    def get_station_distances(self, site_in, useful_sites_in):
        """
        Obtains the distances between stations, determining which are closest to our
        station of interest.
        
        Args:
            site_in               (str):  site_id for our station of interest
            useful_sites_in (list, str):  site_id's for our reference stations
        
        Returns:
            station_distances: pandas.Dataframe containing (sorted) list of distances to stations
                               these will be sorted by distance, with the station of interest removed
                Required Index:
                    site_id   (str): identifiers for the stations
                Required Columns:
                    distance (float): distance to listed station from station of interest, in km

        Assert:
            self.station_data is not None
        """
        if not isinstance(site_in, str):
            raise ValueError('site_in must be a valid string')
        if not isinstance(useful_sites_in, list):
            ValueError('site_in must be a valid list')
        if not all(isinstance(item, str) for item in useful_sites_in):
            raise ValueError('all elements in useful_sites_in list must be strings')

        assert self.station_data is not None, "self.station_data must not be None"

        station_location = (self.station_data.loc[site_in]['latitude'], self.station_data.loc[site_in]['longitude'])
        station_distances = self.calc_station_distances(stations_in=self.station_data.loc[useful_sites_in],
                                                        location=station_location)

        # sort by distance, then drop any station which is the same location as our site of interest
        station_distances = station_distances.sort_values(by='distance', ascending=True)
        station_distances[station_distances.distance==0]=np.nan
        station_distances = station_distances.dropna()

        return station_distances

    def station_listing(self, grouped_data_in):
        """
        Calculates the lists of required sites (those with more data than the minimum required data)
        and reference sites (those with more data than that required for reference purposes).
        
        The requirements for data are defined as the number of days which have at least one reading.
        self.min_years - this is the requirement for the required sites
        self.min_years_reference - this is the requirement for the reference sites
        
        Args:
            grouped_data_in: pandas series object
                Required MultiIndex:
                    site_id: (level 0)
                    timestamp: (level 1)
                Required data:
                    daily count of measurement data (should be in range 0-24 for hourly data)
                    
        Returns:
            required_site_list (list of strings):
                list of sites with a data count > min_years
            useful_site_list (list of strings):
                list of sites with a data count > min_years_reference

        """

        if grouped_data_in is None:
            raise ValueError('grouped_data_in can not be None')
        if not isinstance(grouped_data_in, pd.Series):
            raise ValueError('grouped_data_in must be a pandas series')
        if len(grouped_data_in.index.names) != 2:
            raise ValueError(
            'grouped_data_in series index must have two levels: site_id (level 0) and timestamp (level 1)')


        site_list_interior = grouped_data_in.index.levels[0]
        required_site_list = []
        useful_site_list = []

        for site in site_list_interior:
            try:
                date_num = len(grouped_data_in.loc[(site,),])
            except:
                date_num = 0
            if date_num > self.min_years * 365:
                required_site_list.append(site)
                print('\t{} has {} years of data'.format(site, date_num / 365))
            if date_num > self.min_years_reference * 365:
                useful_site_list.append(site)

        return required_site_list, useful_site_list


    ### function for creating date objects
    def parse_calcs_date(self, date_in):
        return datetime.strptime(date_in, self.DATE_CALCS_FORMAT)

    def calc_nanmean(self, data_in):
        return np.nanmean(data_in)
