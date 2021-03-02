import unittest
from os import path
import pandas as pd
import numpy as np

from environmental_data_modules.met_post_processor import MetPostProcessor


class TestMetPostProcessor(unittest.TestCase):
    """
    Tests for the MetPostProcessor class
    """

    def setUp(self):
        dir, _ = path.split(__file__)
        self.out_dir = path.join(dir, 'output')
        self.verbose = 0
        self.file_in = path.join(dir, 'data', 'processing', 'MEDMI',
                                 'Met_extracted_temp_extras-rel_hum-pressure-dewpoint_reduced.csv')
        self.stations_in_file = path.join(dir, 'data', 'processing', 'stations_in.csv')
        self.result_filename = path.join(dir, 'data', 'OK', 'results_AURN', 'result_calc_station_distance.csv')
        self.stations_in = pd.read_csv(self.stations_in_file)
        self.df_bad_stations_in = pd.DataFrame({'bad_col': [8, 6, 2, 9]})
        self.bad_list_params = [None, 'bad', 20, [10, 'mixed list', 8.3, self.df_bad_stations_in],
                                self.df_bad_stations_in]
        self.result = pd.read_csv(path.join(dir, 'data', 'OK', 'results_MEDMI', 'Met_ppd_daily_mean_max_temp_RH_pres.csv'),
                                  parse_dates=['timestamp'],
                                  index_col=['timestamp', 'site_id'])

    def test_load_good_params(self):
        """
        Test that an PostProcessor can be initialised with OK parameters
        """
        post_processor = MetPostProcessor(out_dir=self.out_dir, station_data_filename=self.stations_in_file,
                                          verbose=self.verbose)
        self.assertIsNotNone(post_processor)
        self.assertIsInstance(post_processor, MetPostProcessor)


    def test_load_bad_data(self):
        """
        Test that a PostProcessor object can not be initialized with bad parameters
        """
        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                post_processor = MetPostProcessor(self.out_dir, self.stations_in_file, self.verbose)
                post_processor = MetPostProcessor(bad_param, self.stations_in_file, self.verbose)
                post_processor = MetPostProcessor(self.out_dir, self.stations_in_file, bad_param)

    def test_process_ok_params(self):
        """
        Test that a MetPostProcessor.process can be called with OK params
        """
        post_processor = MetPostProcessor(out_dir=self.out_dir, station_data_filename=self.stations_in_file,
                                          verbose=self.verbose)
        with self.assertRaises(Warning):
            post_processor.process(in_file=self.file_in,
                                   date_range=["2017-01-01_00", "2017-05-31_23"])

        post_processor.impute_method_setup()
        result = post_processor.process(in_file=self.file_in,
                                        date_range=["2017-01-01_00", "2018-12-31_23"],
                                        impute_data=False)

        # Compare with model result
        pd.testing.assert_frame_equal(result.round(2), self.result.round(2), check_dtype=False, check_frame_type=False, check_exact=False,
                                      check_less_precise=2)

    def test_process_bad_params(self):
        """
        Test that a MetPostProcessor.process cannot be called with bad params
        """
        post_processor = MetPostProcessor(out_dir=self.out_dir, station_data_filename=self.stations_in_file,
                                          verbose=self.verbose)
        post_processor.impute_method_setup()

        with self.assertRaises(ValueError):
            for bad_param in self.bad_list_params:
                post_processor.process(in_file=bad_param,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=bad_param,
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       impute_data=False)
                post_processor.process(in_file=self.file_in,
                                       date_range=["2017-01-01_00", "2018-12-31_23"],
                                       impute_data=bad_param)
