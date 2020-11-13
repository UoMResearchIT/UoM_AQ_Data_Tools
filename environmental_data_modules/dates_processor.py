from abc import ABCMeta, abstractmethod
from datetime import datetime

class DatesProcessor(object):
    """
        Abstract class, parent of date processors (each implementing a different way of representing dates).
        Requires datetime
    """
    __metaclass__ = ABCMeta


class DateRangeProcessor(DatesProcessor):
    INPUT_DATE_FORMAT = '%Y-%m-%d_%H'
    START_DATE = datetime(2016, 1, 1, 0)

    def __init__(self):
        self.__start = None
        self.__end = None

    @property
    def date_range(self):
        return [self.start, self.end]

    @date_range.setter
    def date_range(self, date_range):
        self.set_start_date(date_range[0])
        self.set_end_date(date_range[1])
        if self.start >= self.end:
            raise ValueError('Start date is not earlier than end date.')

    @property
    def start(self):
        return self.__start

    @property
    def end(self):
        return self.__end

    @staticmethod
    def get_available_dates():
        return [DateRangeProcessor.get_available_start(), DateRangeProcessor.get_available_end]

    @staticmethod
    def get_available_start():
        return DateRangeProcessor.START_DATE

    @staticmethod
    def get_available_end():
        cur_year = datetime.now().year
        # Todo: Doug, when is best available end date to set, based on now() ? If it depends on Met or Aurn, we'll need
        #  to sub-type this class
        #  Currently set to end of previous year.
        return datetime(cur_year - 1, 12, 31, 23)

    def set_start_date(self, date_start):
        if date_start is None:
            self.__start = None
            return
        if not isinstance(date_start, datetime):
            raise TypeError('Start date is not in datetime format')
        if date_start < DateRangeProcessor.get_available_start():
            raise ValueError('Start date is less than minimum date available.')
        self.__start = date_start

    def set_end_date(self, date_end):
        if date_end is None:
            self.__end = None
            return
        if not isinstance(date_end, datetime):
            raise TypeError('End date is not in datetime format')
        if date_end > DateRangeProcessor.get_available_end():
            raise ValueError('End date is greater than maximum date available.')
        self.__end = date_end


class DateYearsProcessor(DatesProcessor):

    START_YEAR = 2016

    def __init__(self):
        self.__years = None

    @property
    def years(self):
        return self.__years

    @years.setter
    def years(self, years):
        try:
            years = set(list(years))
        except Exception:
            raise TypeError('years must be a list. Current input: {}'.format(years))

        available_years = DateYearsProcessor.get_available_years()
        error_years = set(years) - set(available_years)
        assert len(error_years) == 0, \
            "Each year must be contained in available years: {}. Error years: {}".format(
                available_years, str(error_years))
        self.__years = years

    @property
    def available_start(self):
        return DateYearsProcessor.get_available_start()

    @property
    def available_end(self):
        return DateYearsProcessor.get_available_end()

    @staticmethod
    def get_available_years():
        return [yr for yr in range(DateYearsProcessor.get_available_start(), DateYearsProcessor.get_available_end()+1)]

    @staticmethod
    def get_available_start():
        return DateYearsProcessor.START_YEAR

    @staticmethod
    def get_available_end():
        cur_year = datetime.now().year
        # Todo: Doug, when is best available end date to set, based on now() ? If it depends on Met or Aurn, we'll need
        #  to sub-type this class
        #  Currently set to previous year.
        return cur_year-1
