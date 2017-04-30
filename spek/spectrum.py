import numpy as np
import astropy.cosmology


class Spectrum():
    """
    Generic spectrum class

    Spectrum class storing supenova spectral data and enabling basic operations
    including redshifting and photometry estimation

    Parameters
    ----------
    wavelength : array-like
        Wavelength of the corresponding flux values. Must be 1d, float and
        same length as flux

    flux : array-like
        Flux for a corresponding wavelength values. Must be 1d, float and
        same length as wavelength.

    redshift : float, optional
        Redshift of the supernova. If no value is provided z=0 will be assumed
        with the luminosity_distance set to 10pc.

    cosmology : `astropy.cosmology.FlatLambdaCDM`, optional
        Cosmology object used to compute luminosity distances. If no value is
        provided Planck15 cosmology will be assumed by default.
    """

    def __init__(self, wavelength, flux, redshift=0, cosmology=None):
        # Must find a working way of checking if data is an array
        try:
            iter(wavelength)
            iter(flux)
            wavelength = np.array(wavelength)
            flux = np.array(flux)
        except TypeError:
            print('wavelength and flux values must be 1d array-like objects')

        try:
            wavelength = wavelength.astype(float)
            flux = flux.astype(float)
        except ValueError:
            print('wavelength and flux arrays must numeric')

        if not wavelength.shape == flux.shape:
            raise ValueError('wavelength and flux must be have same shape')

        if not wavelength.shape[0] == wavelength.size:
            raise ValueError('wavelength and flux must be 1d arrays')

        if cosmology is None:
            self.__cosmology = astropy.cosmology.Planck15
        elif type(cosmology) == astropy.cosmology.core.FlatLambdaCDM:
            self.__cosmology = cosmology
        else:
            raise TypeError("""cosmology must be of type
                            astropy.cosmology.core.FlatLambdaCDM""")

        self._wavelength = wavelength
        self._flux = flux
        self._z = redshift

    def adjust_redshift(self, new_redshift):
        """
        Move the spectrum to a new redshift

        Moving the redshit changes the wavelength and flux by factors of
        (1 + new_redshift) / (1 + old_redshift). Flux is also adjusted to the
        corresponding luminosity distance.

        Parameters
        ----------
        new_redshift: float
            New redshift for the supernova

        Returns
        -------
        None
        """
        redshift_shift = (1 + new_redshift) / (1 + self._z)
        self._wavelength *= redshift_shift
        self._flux /= redshift_shift
        self._z = new_redshift
