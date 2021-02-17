import unittest
from os import path

from environmental_data_modules.extractor import Extractor


class TestExtractor(unittest.TestCase):
  """
  Tests for the Extractor class
  """

  def setUp(self):
    pass


  def test_direct_init(self):
    """
    Test that an Extractor class can not be initialised diretly
    """