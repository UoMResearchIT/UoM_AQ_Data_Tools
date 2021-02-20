import unittest
from os import path
from datetime import datetime

from environmental_data_modules.dates_processor import DatesProcessor, DateRangeProcessor, DateYearsProcessor

class TestDatesProcessor(unittest.TestCase):
    """
    Tests for the DatesProcessor class
    """

    def setUp(self):
        self.years = [2017, 2018, 2019]
        self.date_range = [datetime(2010, 1, 1, 0), datetime(2011, 1, 1, 0)]

    def test_init_date_range(self):
        self.date_range_processor = DateRangeProcessor()
        self.assertIsNotNone(self.date_range_processor)
        self.assertIsInstance(self.date_range_processor, DateRangeProcessor)

    def test_init_date_years(self):
        self.date_years_processor = DateYearsProcessor()
        self.assertIsNotNone(self.date_years_processor)
        self.assertIsInstance(self.date_years_processor, DateYearsProcessor)


    def test_set_date_range_ok_params(self):
        """
        Test constructor with correct inputs
        """
        self.date_range_processor = DateRangeProcessor()
        self.date_range_processor.date_range = self.date_range
        self.assertIsNotNone(self.date_range_processor.date_range)
        self.assertIsInstance(self.date_range_processor.date_range, list)
        self.assertTrue(len(self.date_range_processor.date_range) == 2)
        self.assertTrue(self.date_range_processor.start == self.date_range[0])
        self.assertTrue(self.date_range_processor.end == self.date_range[1])


    def test_set_date_range_bad_params(self):
        """
        Test constructor with invalid inputs
        """
        self.date_range_processor = DateRangeProcessor()
        bad_dates = [['bad1', 'bad2'],
                     [0, 3],
                     [datetime(2011, 1, 12, 0), datetime(2010, 1, 2, 5)],
                     [self, self]
                     ]

        with self.assertRaises(AssertionError):
            for bad_dates in bad_dates:
                self.date_range_processor.date_range = bad_dates

    def test_date_range_get_available_dates(self):
        """
        Test get_available_dates method
        """
        self.date_range_processor = DateRangeProcessor()
        available_dates = self.date_range_processor.get_available_dates()
        self.assertIsNotNone(available_dates)
        self.assertIsInstance(available_dates, list)
        self.assertTrue(len(available_dates) == 2)
        self.assertTrue(available_dates[0] <= available_dates[1])

    def test_set_date_years_ok_params(self):
        """
        Test setting years with correct inputs
        """
        self.date_years_processor = DateYearsProcessor()
        self.date_years_processor.years = self.years
        self.assertIsNotNone(self.date_years_processor.years)
        self.assertIsInstance(self.date_years_processor.years, list)
        self.assertEqual(self.date_years_processor.years, self.years)


    def test_set_date_years_bad_params(self):
        """
        Test setting years with invalid inputs
        """
        self.date_years_processor = DateYearsProcessor()
        bad_dates = [['2010', '2016', '2017'],
                     [0, 3],
                     [datetime(2011, 1, 12, 0), datetime(2010, 1, 2, 5)],
                     [self, self]
                     ]

        with self.assertRaises(AssertionError):
            for bad_dates in bad_dates:
                self.date_years_processor.years = bad_dates


    def test_date_years_get_available_dates(self):
        """
        Test get_available_years method
        """
        self.date_years_processor = DateYearsProcessor()
        available_years = self.date_years_processor.get_available_years()
        self.assertIsNotNone(available_years)
        self.assertIsInstance(available_years, list)
        self.assertTrue(available_years <= sorted(available_years))

