import unittest
from os import path
from datetime import datetime

from environmental_data_modules.regions_processor import RegionsProcessor, RegionRectProcessor, RegionPolyProcessor

class TestRegionsProcessor(unittest.TestCase):
    """
    Tests for the RegionsProcessor class
    """

    def setUp(self):
        self.longitude_range = [-15, 15]
        self.latitude_range = [40, 50]

    def test_init_region_rect_processor(self):
        self.region_rect_processor = RegionRectProcessor()
        self.assertIsNotNone(self.region_rect_processor)
        self.assertIsInstance(self.region_rect_processor, RegionRectProcessor)

    def test_init_region_poly_years(self):
        self.region_poly_processor = RegionPolyProcessor()
        self.assertIsNotNone(self.region_poly_processor)
        self.assertIsInstance(self.region_poly_processor, RegionPolyProcessor)


    def test_set_region_rect_ok_params(self):
        """
        Test constructor with correct inputs
        """
        self.region_rect_processor = RegionRectProcessor()
        self.region_rect_processor.set_region(self.latitude_range, self.longitude_range)
        self.assertIsNotNone(self.region_rect_processor.latitude_range)
        self.assertIsNotNone(self.region_rect_processor.longitude_range)
        self.assertIsInstance(self.region_rect_processor.latitude_range, list)
        self.assertIsInstance(self.region_rect_processor.longitude_range, list)
        self.assertTrue(len(self.region_rect_processor.longitude_range) == 2)
        self.assertTrue(len(self.region_rect_processor.latitude_range) == 2)

    def test_set_region_rect_bad_params(self):
        """
        Test constructor with invalid inputs
        """
        self.region_rect_processor = RegionRectProcessor()

        with self.assertRaises(AssertionError):
            self.region_rect_processor.set_region(['bad', 'bad'], ['bad', 'bad'])
            self.region_rect_processor.set_region([50, 69], [1, 4])
            self.region_rect_processor.set_region([self, self], [self, self])
            self.region_rect_processor.set_region([True, False], [False, True])
