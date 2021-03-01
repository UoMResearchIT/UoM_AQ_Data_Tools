import unittest
from os import path


from environmental_data_modules.post_processor import PostProcessor


class TestPostProcessor(unittest.TestCase):
    """
    Tests for the PostProcessor class
    """

    def setUp(self):
        dir, _ = path.split(__file__)
        self.out_dir = path.join(dir, 'output')
        self.verbose = 0

    def test_direct_init(self):
        """
        Test that a PostProcessor class can not be initialised and then called directly
        """
        with self.assertRaises(NotImplementedError):
            post_processor = PostProcessor(self.out_dir, self.verbose)