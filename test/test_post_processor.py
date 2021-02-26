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
        Test that an PostProcessor class can not be initialised and then called directly
        """
        post_processor = PostProcessor(self.out_dir, self.verbose)
        with self.assertRaises(NotImplementedError):
            post_processor.process(self.file_in)

