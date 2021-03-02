import unittest
from os import path
import pandas as pd
import numpy as np

from environmental_data_modules.aurn_post_processor import AurnPostProcessor


class TestAurnPostProcessor(unittest.TestCase):
    """
    Tests for the AurnPostProcessor class
    """

    def setUp(self):
        dir, _ = path.split(__file__)
        self.metadata_filename = path.join(dir, 'data', 'extraction', 'AURN_metadata.RData')
        self.metadata_url = 'https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData'
        self.out_dir = path.join(dir, 'output')
        self.verbose = 0
        self.file_in = path.join(self.out_dir, 'AURN_extracted_test.csv')
        self.stations_in_file = path.join(dir, 'data', 'processing', 'stations_in.csv')
        self.stat_locn = (51.91614, -0.99958)
        self.result_filename = path.join(dir, 'data', 'OK', 'results_AURN', 'result_calc_station_distance.csv')
        self.site_in = '61'
        self.useful_sites_in = ['50', '80', '500', '800']
        self.stations_in = pd.read_csv(self.stations_in_file)
        self.df_bad_stations_in = pd.DataFrame({'bad_col': [8, 6, 2, 9]})
        self.bad_list_params = [None, 'bad', 20, [10, 'mixed list', 8.3, self.df_bad_stations_in],
                                self.df_bad_stations_in]
        self.site_list = ['CARD', 'LEED', 'ED3', 'NOTT', 'LEAM', 'MY1', 'LON6', 'CHBO']
        self.species_list = ['O3', 'PM10', 'PM2.5', 'NO2', 'NOXasNO2', 'SO2']
        self.result = pd.read_csv(path.join(dir, 'data', 'OK', 'results_AURN', 'aurn_post_process_result.csv'),
                                  parse_dates=['timestamp'],
                                  index_col=['timestamp', 'site_id'])

    def test_load_good_params(self):
        """
        Test that an PostProcessor can be initialised with OK parameters
        """
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)
        self.assertIsNotNone(post_processor)
        self.assertIsInstance(post_processor, AurnPostProcessor)


    def test_load_bad_data(self):
        """
        Test that a PostProcessor object can not be initialized with bad parameters
        """
        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                post_processor = AurnPostProcessor(bad_param, self.out_dir, self.verbose)
                post_processor = AurnPostProcessor(self.metadata_filename, bad_param, self.verbose)
                post_processor = AurnPostProcessor(self.metadata_filename, self.out_dir, bad_param)

    def test_calc_station_distances_ok_params(self):
        """
        Test that a PostProcessor.calc_station_distances
        """

        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)
        result = post_processor.calc_station_distances(self.stations_in, self.stat_locn)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        stored_result = pd.read_csv(self.result_filename, index_col='site_id')
        self.assertTrue(result.round(3).equals(stored_result.round(3)))

    def test_calc_station_distances_bad_params(self):
        """
        Test that a PostProcessor.calc_station_distances fails with bad params
        """
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)

        with self.assertRaises(ValueError):
            for bad_param in self.bad_list_params:
                post_processor.calc_station_distances(stations_in=bad_param,
                                                      location=self.stat_locn)
                post_processor.calc_station_distances(stations_in=self.stations_in,
                                                      location=bad_param)

    def test_get_station_distances_ok_params(self):
        """
        Test that a PostProcessor.get_station_distances can be called with OK params
        """
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)
        with self.assertRaises(AssertionError):
            result = post_processor.get_station_distances(self.site_in, self.useful_sites_in)


    def test_get_station_distances_bad_params(self):
        """
        Test that a PostProcessor.get_station_distances fails with bad params
        """

        bad_list_params = [None, 'bad', 20, [10, 'mixed list', 8.3, self]]  # Bad lists
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)

        with self.assertRaises(ValueError):
            for bad_param in bad_list_params:
                post_processor.get_station_distances(site_in=bad_param,
                                                     useful_sites_in=self.useful_sites_in)
                post_processor.get_station_distances(site_in=self.site_in,
                                                     useful_sites_in=bad_param)

    def test_station_listing_ok_params(self):
        """
        Test that a PostProcessor.listing can be called with OK params
        """
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)

        idx = pd.MultiIndex.from_product([['1', '200', '400', '420', '555'], ['2017-04-01', '2017-04-02']])
        data = np.random.randn(10).round()
        series = pd.Series(data=data, index=idx)

        required, useful = post_processor.station_listing(series)
        self.assertIsNotNone(required)
        self.assertIsNotNone(useful)
        self.assertIsInstance(required, list)
        self.assertIsInstance(useful, list)

    def test_station_listing_bad_params(self):
        """
        Test that a PostProcessor.station_listing fails with bad params
        """
        idx = pd.MultiIndex.from_product([['1', '200', '400', '420', '555'], ['2017-04-01', '2017-04-02']])
        idx_1 = pd.MultiIndex.from_product([[1, 3, 6, 9, 10], ['2017-04-01', '2017-04-02']])
        idx_2 = pd.MultiIndex.from_product([['1', '200', '400', '420', '555'], [10, 848.98]])
        data = np.random.randn(10).round()
        data_2 = np.random.randn(10)
        bad_series_1 = pd.Series(data=data, index=idx_1)
        bad_series_2 = pd.Series(data=data, index=idx_2)
        bad_series_3 = pd.Series(data=data_2, index=idx)

        bad_list_params = [None, bad_series_1, bad_series_2, bad_series_3, 'bad', 20, [10, 'mixed list', 8.3, self]]
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)

        with self.assertRaises(ValueError):
            for bad_param in bad_list_params:
                post_processor.station_listing(bad_param)

    def test_process_ok_params(self):
        """
        Test that a AurnPostProcessor.process can be called with OK params
        """
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)

        with self.assertRaises(Warning):
            post_processor.process(in_file=self.file_in,
                                   date_range=["2017-01-01_00", "2017-05-31_23"],
                                   site_list=self.site_list,
                                   species_list=['O3'])

        post_processor.impute_method_setup()
        result = post_processor.process(in_file=self.file_in,
                                        date_range=["2017-01-01_00", "2018-12-31_23"],
                                        site_list=self.site_list,
                                        species_list=self.species_list,
                                        impute_data=False)

        # Compare with model result
        pd.testing.assert_frame_equal(result, self.result, check_dtype=False, check_frame_type=False, check_exact=False)

    def test_process_bad_params(self):
        """
        Test that a AurnPostProcessor.process cannot be called with bad params
        """
        post_processor = AurnPostProcessor(self.metadata_filename, out_dir=self.out_dir, verbose=self.verbose)
        post_processor.impute_method_setup()

        with self.assertRaises(ValueError):
            for bad_param in self.bad_list_params:
                post_processor.process(in_file=bad_param,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       site_list=self.site_list,
                                       species_list=self.species_list,
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=bad_param,
                                       site_list=self.site_list,
                                       species_list=self.species_list,
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       site_list=bad_param,
                                       species_list=self.species_list,
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       site_list=self.site_list,
                                       species_list=bad_param,
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       site_list=self.site_list,
                                       species_list=self.species_list,
                                       impute_data=bad_param)
