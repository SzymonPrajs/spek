import os
import numpy as np


class Filters():
    """
    Filter responces class

    Generic functionality for loading, providing statistics and using
    astronomical filter responces. This can be used completely independantly of
    of other classes that are part of this package and should work with both
    numpy array and pandas DataFrames

    Parameters
    ----------
    filters_directory : str
        Path to a directory that stores filter responces. The files must have
        a *.dat extention and contain two culumns: wavelength and bandpass.

    load_all : bool, optional
        Filter respoce files are not loaded by default to minimise the load
        time for the module in case a large number of filter responces are
        provided. This can be changed by setting this flag to True.
    """
    def __init__(self, filters_directory, load_all=False):
        if os.path.exists(filters_directory):
            self._filter_directory = filters_directory
        else:
            raise IOError('File {} does not exist!'.format(filters_directory))

    def make_filter_list(self):
        """
        Scan though the filter_drectory and sets the internal list of filter
        responces that can be later loaded using the load_filters method.
        """
        pass

    def load_filters(self, filter_name=None, load_all=False):
        """
        """
        # A lot to do here still
        np.loadtxt(os.path.join(self._filters_directory,
                                filter_name), unpack=True)

    def list_available(self):
        """
        List all available filters
        """
        pass

    def list_loaded(self):
        """
        List all filters that have already been loaded
        """
        pass
