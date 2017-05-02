import os
import glob
import numpy as np
from astropy.constants import c


class _Filter():
    """
    Filter response data structure (for internal use only)

    Stores individual filter responces and their statistics. The objects are
    loaded and accessed by the user though the main Filters class.

    Parameters
    ----------
    filter_path : str
        Path to the filter response to be loaded
    """
    def __init__(self, filter_path):
        # TODO: Do all appropriate checks
        self.file_path = filter_path
        self.base_name = os.path.basename(self.file_path)
        self.name = os.path.splitext(self.base_name)[0]

        # TODO: Check if file has correct format
        filter_data = np.loadtxt(filter_path, unpack=True)
        self.wavelength = filter_data[0]
        self.bandpass = filter_data[1]
        self.area = np.trapz(self.bandpass, x=self.wavelength)
        jansky = 3631 * 1e-23 * c.value * 1e10 / (self.wavelength)**2
        flux = jansky * self.bandpass
        self.ab_zero_point = -2.5*np.log10(np.trapz(flux, x=self.wavelength) /
                                           self.area)
        self.central_wavelength = np.trapz(self.bandpass * self.wavelength,
                                           x=self.wavelength) / self.area


class Filters():
    """
    Filter responces utility

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
            raise IOError('Path does not exist: {}'.format(filters_directory))

        self.make_filter_list()
        self.filters = {}
        self._loaded_filters = []

        if load_all is True:
            self.load_filters(load_all=True)

    def __getitem__(self, filter_name):
        """
        `Official` access point for indivisual filter

        For ease and readibility the __getitem__ method was overwritten to
        returns an object of the class _Filter after checking if a given
        filter has been loaded. This is safer than the accessing the
        Filters.filters dictionary dirrectly however it is marginally slower.

        Parameters
        ----------
        filter_name : str
            Name of the filter to be accessed

        Returns
        -------
        filter : _Filter
            Object of the class _Filter containing filter responces and stats

        Example
        -------
        >>> f = Filters('path/to/filters', load_all=True)
        >>> f.filters['DES_g'].name
        'DES_g'
        >>> f['DES_g'].name
        'DES_g'
        """
        return self.filters[filter_name]

    def make_filter_list(self):
        """
        Scan though the filter_drectory and sets the internal list of filter
        responces that can be later loaded using the load_filters method.
        """
        self.__path_list = np.array(glob.glob(self._filter_directory+'/*.dat'))
        self.__base_names = list(map(os.path.basename, self.__path_list))
        self.__filter_names = np.array([b[0] for b in map(os.path.splitext,
                                                          self.__base_names)])

    def load_filters(self, filter_name=None, load_all=False):
        """
        """
        if load_all is True:
            filters_to_load = np.array(self.__filter_names)
        elif filter_name is not None:
            filters_to_load = np.array(filter_name)
        else:
            raise ValueError('List of filters was not provided')

        for ftl in filters_to_load:
            idx = np.where(self.__filter_names == ftl)
            if idx[0].size == 0:
                raise Warning('No filter responce for {}'.format(ftl))
            else:
                self.filters[ftl] = _Filter(self.__path_list[idx][0])
                self._loaded_filters.append(ftl)

    def list_available(self):
        """
        List all available filters
        """
        return self.__filter_names

    def list_loaded(self):
        """
        List all filters that have already been loaded
        """
        return np.unique(np.array(self._loaded_filters))
