from abc import ABCMeta, abstractmethod

class RegionsProcessor(object):
    """
        Abstract class, parent of regions processors (each implementing a different way of representing regions).
        Requires datetime
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def set_region(self):
        raise NotImplementedError("Must override set_region")


class RegionRectProcessor(RegionsProcessor):
    """
        Class used for processing region of interest as rectangles, via latitude/longitude ranges
        # Note: If there is another useful way to represent rectangles, (e.g. other co-ordinate system), you would
                make a sub-class of this class and place the current content into (e.g.) RegionLatLongProcessor,
                and the other way represented in (e.g.) Region[Other-way's-name]Processor .
    """

    # Define 'absolute' constants
    UK_LATITUDES = [48., 60.]
    UK_LONGITUDES = [-11., 3.]
    LONGITUDE_RANGE = [-180., 360.]
    LATITUDE_RANGE = [-90., 90.]

    def __init__(self):
        """ Initialise instance of the RegionRectProcessor class.
            Initialises the private class variables

            Returns:
                Initialised instance of RegionRectProcessor

        """
        self.latitude_range = RegionRectProcessor.UK_LATITUDES
        self.longitude_range = RegionRectProcessor.UK_LONGITUDES

    @property
    def latitude_range(self):
        return self.__latitude_range

    @latitude_range.setter
    def latitude_range(self, range):
        """ Sets the latitude range of interest

            Args:
                range: (list of 2 floats) The start (range[0]) and end (range[1]) of the latitude range of interest
                The 2 values must both fall within the allowed global range of latitude values

            Returns:
                None

        """
        assert isinstance(range, list), 'renge must be a valid list'
        assert len(range) == 2, "range must be a list of 2 values"
        assert isinstance(range[0], float) or isinstance(range[0], int), 'range[0] value must be numeric'
        assert isinstance(range[1], float) or isinstance(range[1], int), 'range[1] value must be numeric'
        val_1 = float(range[0])
        val_2 = float(range[1])
        if not ((RegionRectProcessor.LATITUDE_RANGE[0] <= val_1) and (val_1 <= RegionRectProcessor.LATITUDE_RANGE[1])):
            raise ValueError('Latitude first value falls outside global range')
        if not ((RegionRectProcessor.LATITUDE_RANGE[0] <= val_2) and (val_2 <= RegionRectProcessor.LATITUDE_RANGE[1])):
            raise ValueError('Latitude last value falls outside global range')
        self.__latitude_range = [val_1, val_2]

    @property
    def longitude_range(self):
        return self.__longitude_range

    @longitude_range.setter
    def longitude_range(self, range):
        """ Sets the longitude range of interest

            Args:
                range: (list of 2 floats) The start (range[0]) and end (range[1]) of the longitude range of interest
                The 2 values must both fall within the allowed global range of longitude values

            Returns:
                None

        """
        assert isinstance(range, list), 'renge must be a valid list'
        assert len(range) == 2, "range must be a list of 2 values"
        assert isinstance(range[0], float) or isinstance(range[0], int), 'range[0] value must be numeric'
        assert isinstance(range[1], float) or isinstance(range[1], int), 'range[1] value must be numeric'
        val_1 = float(range[0])
        val_2 = float(range[1])
        if not (RegionRectProcessor.LONGITUDE_RANGE[0] <= val_1 <= RegionRectProcessor.LONGITUDE_RANGE[1]):
            raise ValueError('Longitude first value falls outside global range')
        if not (RegionRectProcessor.LONGITUDE_RANGE[0] <= val_2 <= RegionRectProcessor.LONGITUDE_RANGE[1]):
            raise ValueError('Longitude last value falls outside global range')
        self.__longitude_range = [val_1, val_2]

    def set_region(self, latitude_range, longitude_range):
        """ Sets the latitude range of interest

            Args:
                latitude_range: (list of 2 floats) The start (range[0]) and end (range[1]) of the
                            latitude range of interest
                longitude_range: (list of 2 floats) The start (range[0]) and end (range[1]) of the
                            longitude range of interest
            Returns:
                None

        """
        #Todo: Doug do we need to test if start < end for both long and latitude?
        self.latitude_range = latitude_range
        self.longitude_range = longitude_range


class RegionPolyProcessor(RegionsProcessor):
    """
        Class used for processing list of polygons as region of interest
        Note: Added as a placeholder and example of a RegionsProcessor sub-class
    """
    # Todo

    POLYGONS_UK = ['Fake uk polygons']

    def __init__(self):
        self.polygons = RegionPolyProcessor.POLYGONS_UK

    @property
    def polygons(self):
        return self.__polygons

    @polygons.setter
    def polygons(self, polygons):
        # Test the polygons argument
        self.__polygons = polygons

    def set_region(self, polygon_list):
        # If OK, use the polygon
        self.polygons = polygon_list

