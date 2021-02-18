import unittest
from os import path

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
    self.site_list = []
    self.years = [2017]


  def test_init_ok_params(self):
    """
    Test constructor with correct inputs
    """

    # No defaults, both filename and url
    extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                              metadata_url=self.metadata_url,
                              out_dir=self.out_dir,
                              verbose=0)
    self.assertIsNotNone(extractor)
    self.assertIsInstance(extractor, AurnExtractor)

    # No defaults, filename only
    extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                              metadata_url=None,
                              out_dir=self.out_dir,
                              verbose=0)
    self.assertIsNotNone(extractor)
    self.assertIsInstance(extractor, AurnExtractor)

    # No defaults, URL only
    extractor = AurnExtractor(metadata_filename=None,
                              metadata_url=self.metadata_url,
                              out_dir=self.out_dir,
                              verbose=0)
    self.assertIsNotNone(extractor)
    self.assertIsInstance(extractor, AurnExtractor)

    # Will use default filename (None)
    extractor = AurnExtractor(metadata_url=self.metadata_url,
                              out_dir=self.out_dir,
                              verbose=0)
    self.assertIsNotNone(extractor)
    self.assertIsInstance(extractor, AurnExtractor)

    # Will use default filename (None) and default url
    extractor = AurnExtractor(out_dir=self.out_dir,
                              verbose=0)
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
                                  verbose=0)

        extractor = AurnExtractor(metadata_url='bad url',
                                  out_dir=self.out_dir,
                                  verbose=0)

        extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                                  out_dir=108748485,
                                  verbose=0)

        extractor = AurnExtractor(metadata_url=self.metadata_filename,
                                  out_dir=self.out_dir,
                                  verbose='bad verbose')


  def test_extract_data_OK_params(self):
    """
    Test extract_data with OK inputs
    """
    extractor = AurnExtractor(metadata_filename=self.metadata_filename,
                              out_dir=self.out_dir,
                              verbose=0)
    extractor.extract_data(
      years=self.years,
      site_list=self.site_list,
      save_to_csv=False,
      outfile_suffix=self.outfile_suffix)


'''
  def test_extract_data_bad_params(self):
    """
    Test extract_data with bad inputs
    """
    extractor = AurnExtractor(metadata_filename=self.load_metadata_filename,
                              metadata_url=self.metadata_url,
                              out_dir=self.out_dir, verbose=0)

    extractor.extract_data(
      years=100000,
      site_list=self.site_list,
      save_to_csv=False,
      outfile_suffix=self.outfile_suffix)

    extractor.extract_data(
      years=[2017],
      site_list=10,
      save_to_csv=False,
      outfile_suffix=self.outfile_suffix)

    extractor.extract_data(
      years=2017,
      site_list=self.site_list,
      save_to_csv=False,
      outfile_suffix=10)'''


