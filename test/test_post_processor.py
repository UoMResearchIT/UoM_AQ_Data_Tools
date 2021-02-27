import unittest
from os import path
import pandas as pd

from environmental_data_modules.post_processor import PostProcessor


class TestPostProcessor(unittest.TestCase):
    """
    Tests for the PostProcessor class
    """

    def setUp(self):
        dir, _ = path.split(__file__)
        self.out_dir = path.join(dir, 'output')
        self.verbose = 0
        self.file_in = path.join(self.out_dir, 'AURN_extracted_test.csv')
        self.stations_in_file = path.join(dir, 'data', 'processing', 'stations_in.csv')
        self.stat_locn = (51.91614, -0.99958)
        self.result_filename = path.join(dir, 'data', 'OK', 'results_AURN', 'result_calc_station_distance.csv')
        self.site_in = 61
        self.useful_sites_in = [50, 80, 500, 800]

    def test_load_good_params(self):
        """
        Test that an PostProcessor can be initialised with OK parameters
        """
        post_processor = PostProcessor(self.out_dir, self.verbose)
        self.assertIsNotNone(post_processor)
        self.assertIsInstance(post_processor, PostProcessor)


    def test_load_bad_data(self):
        """
        Test that a PostProcessor object can not be initialized with bad parameters
        """
        with self.assertRaises(AssertionError):
            post_processor = PostProcessor(10, self.verbose)
            post_processor = PostProcessor('?>>... **', self.verbose)
            post_processor = PostProcessor(self.out_dir, 'bad verbose')
            post_processor = PostProcessor(self.out_dir, 80.4)

    def test_direct_init(self):
        """
        Test that a PostProcessor class can not be initialised and then called directly
        """
        post_processor = PostProcessor(self.out_dir, self.verbose)
        with self.assertRaises(NotImplementedError):
            post_processor.process(self.file_in)

    def test_calc_station_distances(self):
        """
        Test that a PostProcessor.calc_station_distances
        """
        stations_in = pd.read_csv(self.stations_in_file)
        post_processor = PostProcessor(self.out_dir, self.verbose)
        result = post_processor.calc_station_distances(stations_in, self.stat_locn)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)
        stored_result = pd.read_csv(self.result_filename, index_col='site_id')
        self.assertTrue(result.round(3).equals(stored_result.round(3)))



    def test_get_station_distances_ok_params(self):
        """
        Test that a PostProcessor.calc_station_distances can be called with OK params
        """
        post_processor = PostProcessor(self.out_dir, self.verbose)
        with self.assertRaises(AssertionError):
            result = post_processor.get_station_distances(self.site_in, self.useful_sites_in)

    def test_get_station_distances_bad_params(self):
        """
        Test that a PostProcessor.calc_station_distances fails with bad params
        """

        bad_list_params = ['bad', 20, [10, 'mixed list', 8.3, self]]  # Bad lists
        post_processor = PostProcessor(self.out_dir, self.verbose)

        with self.assertRaises(ValueError):
            for bad_param in bad_list_params:
                post_processor.get_station_distances(site_in=bad_param,
                                                     useful_sites_in=self.useful_sites_in)
                post_processor.get_station_distances(site_in=self.site_in,
                                                     useful_sites_in=bad_param)