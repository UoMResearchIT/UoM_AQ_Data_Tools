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
            raise Warning('Imputation requested but imputer is currently None')
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

    def calc_station_distances(self, stations_in, stat_location):
        """
        Calculates the distances between stations.
        
        Args:
            stations_in: pandas.Dataframe containing station locations
                Required Index:
                    SiteID   (str): identifiers for the stations
                Required Columns:
                    Latitude (float):
                    Longitude (float):
            stat_location (tuple, float): (latitude, longitude) of station of interest
        
        Returns:
            station_distances: pandas.Dataframe containing (sorted) list of distances to stations
                Required Index:
                    SiteID   (str): identifiers for the stations
                Required Columns:
                    Distance (float): distance to listed station from station of interest, in km 
        """
        station_distances = pd.DataFrame(index=stations_in.index)
        station_distances['Distance'] = np.nan

        for index, row in stations_in.iterrows():
            new_location = (row['Latitude'],row['Longitude'])

            station_distances.loc[index]['Distance'] = distance.distance(stat_location,new_location).km

        return station_distances


    def get_station_distances(self, site_in, useful_sites_in):
        """
        Obtains the distances between stations, determining which are closest to our
        station of interest.
        
        Args:
            site_in               (str):  siteID for our station of interest
            useful_sites_in (list, str):  siteID's for our reference stations
        
        Returns:
            station_distances: pandas.Dataframe containing (sorted) list of distances to stations
                               these will be sorted by distance, with the station of interest removed
                Required Index:
                    SiteID   (str): identifiers for the stations
                Required Columns:
                    Distance (float): distance to listed station from station of interest, in km 
        """

        station_location = (self.station_data.loc[site_in]['Latitude'], self.station_data.loc[site_in]['Longitude'])
        station_distances = self.calc_station_distances(stations_in=self.station_data.loc[useful_sites_in], \
                                                   stat_location=station_location)

        # sort by distance, then drop any station which is the same location as our site of interest
        station_distances = station_distances.sort_values(by='Distance',ascending=True)
        station_distances[station_distances.Distance==0]=np.nan
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
                    SiteID: (level 0)
                    Date: (level 1)
                Required data:
                    daily count of measurement data (should be in range 0-24 for hourly data)
                    
        Returns:
            required_site_list (list of strings):
                list of sites with a data count > min_years
            useful_site_list (list of strings):
                list of sites with a data count > min_years_reference
        """

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
