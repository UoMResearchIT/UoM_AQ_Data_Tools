# -*- coding: utf-8 -*-

import unittest
from os import path
import pandas as pd
from pandas.util.testing import assert_frame_equal
from datetime import datetime

from environmental_data_modules.met_extractor import MetExtractor


class TestMetExtractor(unittest.TestCase):
    """
    Tests for the Extractor class
    """

    def setUp(self):
        dir, _ = path.split(__file__)
        self.out_dir = path.join(dir, 'output')
        self.verbose = 0
        self.measurements = ['temperature', 'rain', 'wind', 'alnus', 'urtica']
        self.result = {'alnus': pd.read_csv(path.join(dir, 'data', 'OK', 'results_MEDMI',
                                                      'Met_extracted_pollen-alnus_reduced.csv'),
                                            parse_dates=['timestamp'], index_col='timestamp', skiprows=1),
                       'urtica': pd.read_csv(path.join(dir, 'data', 'OK', 'results_MEDMI',
                                                       'Met_extracted_pollen-urtica_reduced.csv'),
                                             parse_dates=['timestamp'], index_col='timestamp', skiprows=1),
                       'rain': pd.read_csv(path.join(dir, 'data', 'OK', 'results_MEDMI',
                                                     'Met_extracted_rain_reduced.csv'),
                                           parse_dates=['timestamp'], index_col='timestamp', skiprows=1),
                       'temp': pd.read_csv(path.join(dir, 'data', 'OK', 'results_MEDMI',
                                                     'Met_extracted_temp_extras-rel_hum-pressure-dewpoint_reduced.csv'),
                                           parse_dates=['timestamp'], index_col='timestamp', skiprows=1),
                       'wind': pd.read_csv(path.join(dir, 'data', 'OK', 'results_MEDMI',
                                                     'Met_extracted_wind_reduced.csv'),
                                           parse_dates=['timestamp'], index_col='timestamp', skiprows=1)}

        self.bad_list_params = [[], 'bad', 20, [10, 'mixed list', 8.3, self]]  # Bad lists

        self.extractor = MetExtractor(out_dir=self.out_dir,
                                      verbose=0)

        self.latitudes = [53, 55]
        self.longitudes = [-5, -3]
        self.date_range = ['2017-06-01_0', '2017-06-30_23']
        self.outfile_suffix = "reduced"
        self.extra_measurements = True

    def test_init_ok_params(self):
        """
        Test constructor with correct inputs
        """

        # No defaults, both filename and url
        extractor = MetExtractor(out_dir=self.out_dir,
                                 verbose=self.verbose)
        self.assertIsNotNone(extractor)
        self.assertIsInstance(extractor, MetExtractor)


    def test_init_bad_params(self):
        """
        Test constructor with bad inputs
        """
        with self.assertRaises(AssertionError):
            extractor = MetExtractor(out_dir=108748485,
                                     verbose=self.verbose)
            extractor = MetExtractor(out_dir='&??. ..£',
                                     verbose=self.verbose)
            extractor = MetExtractor(out_dir=[self.out_dir],
                                     verbose=self.verbose)

            extractor = MetExtractor(out_dir=self.out_dir,
                                     verbose='bad verbose')
            extractor = MetExtractor(out_dir=self.out_dir,
                                     verbose=[3])


    def test_extract_data_OK_params(self):
        """
        Test extract_data with OK inputs
        """
        if MetExtractor.dataset is not None and callable(MetExtractor.dataset):
            for measurement in self.measurements:
                class_ = MetExtractor.get_class_from_measurement_name(measurement)
                met_extractor = class_(self.out_dir, self.verbose)
                result = met_extractor.extract_data(
                    extract_extra_datasets=self.extra_measurements,
                    latitude_range=self.latitudes,
                    longitude_range=self.longitudes,
                    outfile_suffix=self.outfile_suffix,
                    date_range=self.date_range)

                self.assertIsNotNone(result)
                self.assertIsInstance(result, pd.DataFrame)

                # Compare with model result
                result.set_index('timestamp', inplace=True)
                self.assertTrue(result.round(3).equals(self.result[measurement].round(3)))
        else:
            print('\nCan\'t see MEDMI Dataset object, so could not run MEDMI data extraction tests.')

    def test_extract_data_bad_extra_measurements(self):
        """
        Test extract_data with bad extra measurements param
        """
        if not self.extractor:
            self.extractor = MetExtractor(out_dir=self.out_dir,
                                          verbose=self.verbose)

        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                    extract_extra_datasets=bad_param,
                    latitude_range=self.latitudes,
                    longitude_range=self.longitudes,
                    outfile_suffix=self.outfile_suffix,
                    date_range=self.date_range)

    def test_extract_data_bad_latitudes_longitudes(self):
        """
        Test extract_data with bad latitude or longitude range param
        """
        if not self.extractor:
            self.extractor = MetExtractor(out_dir=self.out_dir,
                                          verbose=self.verbose)

        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                    extract_extra_datasets=self.extra_measurements,
                    latitude_range=bad_param,
                    longitude_range=self.longitudes,
                    outfile_suffix=self.outfile_suffix,
                    date_range=self.date_range)
                self.extractor.extract_data(
                    extract_extra_datasets=self.extra_measurements,
                    latitude_range=self.latitudes,
                    longitude_range=bad_param,
                    outfile_suffix=self.outfile_suffix,
                    date_range=self.date_range)

    def test_extract_data_bad_outfile_suffix(self):
        """
        Test extract_data with bad outfile suffix param
        """
        if not self.extractor:
            self.extractor = MetExtractor(out_dir=self.out_dir,
                                          verbose=self.verbose)

        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                    extract_extra_datasets=self.extra_measurements,
                    latitude_range=self.latitudes,
                    longitude_range=self.longitudes,
                    outfile_suffix=bad_param,
                    date_range=self.date_range)

    def test_extract_data_bad_date_range(self):
        """
        Test extract_data with bad date range param
        """
        if not self.extractor:
            self.extractor = MetExtractor(out_dir=self.out_dir,
                                          verbose=self.verbose)

        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                    extract_extra_datasets=self.extra_measurements,
                    latitude_range=self.latitudes,
                    longitude_range=self.longitudes,
                    outfile_suffix=self.outfile_suffix,
                    date_range=bad_param)