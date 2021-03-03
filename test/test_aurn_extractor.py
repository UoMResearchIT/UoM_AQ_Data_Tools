import unittest
from os import path
import pandas as pd

from environmental_data_modules.aurn_extractor import AurnExtractor


class TestAurnExtractor(unittest.TestCase):
    """
    Tests for the Extractor class
    """

    def setUp(self):
        dir, _ = path.split(__file__)
        self.metadata_filename = path.join(dir, 'data', 'extraction', 'AURN_metadata.RData')
        self.metadata_url = 'https://uk-air.defra.gov.uk/openair/R_data/AURN_metadata.RData'
        self.out_dir = path.join(dir, 'output')
        self.outfile_suffix = 'test'
        self.site_list = ['CARD', 'LEED', 'ED3', 'NOTT', 'LEAM', 'MY1', 'LON6', 'CHBO']
        self.years = [2017, 2018]
        self.verbose = 0
        self.result = pd.read_csv(path.join(dir, 'data', 'OK', 'results_AURN', 'aurn_extracted_result.csv'),
                                  parse_dates=['timestamp'], index_col=['timestamp', 'site_id'])
        self.bad_list_params = [[], 'bad', 20, [10, 'mixed list', 8.3, self]]  # Bad lists

        self.extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                                       out_dir=self.out_dir,
                                       verbose=0)

    def test_init_ok_params(self):
        """
        Test constructor with correct inputs
        """

        # No defaults, both filename and url
        extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                                  metadata_url=self.metadata_url,
                                  out_dir=self.out_dir,
                                  verbose=self.verbose)
        self.assertIsNotNone(extractor)
        self.assertIsInstance(extractor, AurnExtractor)

        # No defaults, filename only
        extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                                  metadata_url=None,
                                  out_dir=self.out_dir,
                                  verbose=self.verbose)
        self.assertIsNotNone(extractor)
        self.assertIsInstance(extractor, AurnExtractor)

        # No defaults, URL only
        extractor = AurnExtractor(metadata_filename=None,
                                  metadata_url=self.metadata_url,
                                  out_dir=self.out_dir,
                                  verbose=self.verbose)
        self.assertIsNotNone(extractor)
        self.assertIsInstance(extractor, AurnExtractor)

        # Will use default filename (None)
        extractor = AurnExtractor(metadata_url=self.metadata_url,
                                  out_dir=self.out_dir,
                                  verbose=self.verbose)
        self.assertIsNotNone(extractor)
        self.assertIsInstance(extractor, AurnExtractor)

        # Will use default filename (None) and default url
        extractor = AurnExtractor(out_dir=self.out_dir,
                                  verbose=self.verbose)
        self.assertIsNotNone(extractor)
        self.assertIsInstance(extractor, AurnExtractor)


    def test_init_bad_params(self):
        """
        Test constructor with bad inputs
        """
        with self.assertRaises(AssertionError):
            extractor = AurnExtractor(metadata_filename='bad_filename',
                                      metadata_url=self.metadata_url,
                                      out_dir=self.out_dir,
                                      verbose=self.verbose)

            extractor = AurnExtractor(metadata_url='bad url',
                                      out_dir=self.out_dir,
                                      verbose=self.verbose)

            extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                                      out_dir=108748485,
                                      verbose=self.verbose)

            extractor = AurnExtractor(metadata_url=self.metadata_filename,
                                      out_dir=self.out_dir,
                                      verbose='bad verbose')


    def test_extract_data_OK_params_local_metadata(self):
        """
        Test extract_data with OK inputs
        """
        result = self.extractor.extract_data(
            years=self.years,
            site_list=self.site_list,
            save_to_csv=True,
            outfile_suffix=self.outfile_suffix)

        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)

        # Compare with model result
        result.sort_values(by=['timestamp', 'site_id'], inplace=True)
        result.set_index(['timestamp', 'site_id'], inplace=True)
        self.assertTrue(result.round(3).equals(self.result.round(3)))

    def test_extract_data_OK_params_url_metadata(self):
        """
        Test extract_data with OK inputs
        """
        self.extractor = AurnExtractor(out_dir=self.out_dir,
                                       verbose=self.verbose)

        result = self.extractor.extract_data(
            years=self.years,
            site_list=self.site_list,
            save_to_csv=False,
            outfile_suffix=self.outfile_suffix)

        self.assertIsNotNone(result)
        self.assertIsInstance(result, pd.DataFrame)

        # Compare with model result
        result.sort_values(by=['timestamp', 'site_id'], inplace=True)
        result.set_index(['timestamp', 'site_id'], inplace=True)
        self.assertTrue(result.round(3).equals(self.result.round(3)))

    def test_extract_data_bad_years(self):
        """
        Test extract_data with bad years param
        """
        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                  years=bad_param,
                  site_list=self.site_list,
                  save_to_csv=False,
                  outfile_suffix=self.outfile_suffix)

    def test_extract_data_bad_site_list(self):
        """
        Test extract_data with bad site_list param
        """
        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                  years=self.years,
                  site_list=bad_param,
                  save_to_csv=False,
                  outfile_suffix=self.outfile_suffix)

    def test_extract_data_bad_species_list(self):
        """
        Test extract_data with bad species_list param
        """
        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                    years=self.years,
                    site_list=self.site_list,
                    species_list=bad_param,
                    save_to_csv=False,
                    outfile_suffix=self.outfile_suffix)

    def test_extract_data_bad_save_to_file(self):
        """
        Test extract_data with bad species_list param
        """
        with self.assertRaises(AssertionError):
            for bad_param in self.bad_list_params:
                self.extractor.extract_data(
                    years=self.years,
                    site_list=self.site_list,
                    save_to_csv=bad_param,
                    outfile_suffix=self.outfile_suffix)

    def test_extract_data_bad_outfile_suffix(self):
        """
        Test extract_data with bad outfile_suffix param
        """
        with self.assertRaises(AssertionError):
            self.extractor.extract_data(
              years=self.years,
              site_list=self.site_list,
              save_to_csv=False,
              outfile_suffix=10)

            # Test bad save_to_csv
            self.extractor.extract_data(
              years=self.years,
              site_list=self.site_list,
              save_to_csv=False,
              outfile_suffix='% %..()')
