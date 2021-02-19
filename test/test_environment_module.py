import unittest
from os import path

from environmental_data_modules.environment_module import EnvironmentModule


class TestEnvironmentModule(unittest.TestCase):
  """
  Tests for the EnvironmentModule class
  """

  def setUp(self):
    pass


  def test_direct_init(self):
    """
    Test that an EnvironmentModule can not be initialised diretly
    """