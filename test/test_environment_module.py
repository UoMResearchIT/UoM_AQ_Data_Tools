import unittest
from os import path

from environmental_data_modules.environment_module import EnvironmentModule


class TestEnvironmentModule(unittest.TestCase):
    """
    Tests for the EnvironmentModule class
    """

    def setUp(self):
        dir, _ = path.split(__file__)
        self.out_dir = path.join(dir, 'output')
        self.verbose = 0

    def test_load_good_params(self):
        """
        Test that an EnvironmentModule can not be initialised directly and then be used (calling methods)
        """
        env_mod = EnvironmentModule(self.out_dir, self.verbose)
        self.assertIsNotNone(env_mod)
        self.assertIsInstance(env_mod, EnvironmentModule)


    def test_load_bad_data(self):
        """
        Test that a RegionEstimator object can be initialized with good data.
        Also check that various other initializations happen within the object.
        """
        with self.assertRaises(AssertionError):
            env_mod = EnvironmentModule(10, self.verbose)
            env_mod = EnvironmentModule('?>>... **', self.verbose)
            env_mod = EnvironmentModule(self.out_dir, 'bad verbose')
            env_mod = EnvironmentModule(self.out_dir, 80.4)

