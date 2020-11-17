#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from abc import ABCMeta, abstractmethod

from environmental_data_modules import EnvironmentModule

class Extractor(EnvironmentModule):
    """
        Abstract class, parent of classes used for extracting data from AURN/MEDMI servers.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def extract_data(self, date_range, save_to_file, outfile_suffix):
        raise NotImplementedError("Must override extract_data")
