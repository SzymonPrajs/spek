import numpy as np
from spek.spectrum import Spectrum
import astropy.units as u
import astropy.constants as const

class BlackbodySED(Spectrum):
    """
    Blackbody SED model

    Child class of the Spectrum class populates the wavelength and flux with
    values calculated based on the

    Parameters
    ----------
    temperature : float
        Photospheric temperature of the object. Value must be given in Kelvins.

    radius : float or `astropy.units.quantity.Quantity`
        Radius of the photosphere. If float is provided it will be assumed to
        be in units of cm. Alternatively an astropy.units.quantity.Quantity
        can be passed which will be internally converted into the correct units

    redshift : float, optional
        Redshift at which to place the object. Default assumes restframe (z=0)

    wavelength : array-like, optional
        Wavelength at which the Planck function is to be evaluated.
        If no input is provided, the default will be set to a range of
        1000A - 10000A
    """
    def __init__(self, tempeture, radius, redshift=0, wavelength=None):
        try:
            float(tempeture)
        except:
            raise ValueError('Temperature must be a numeric value')

        if type(radius) != u.quantity.Quantity:
            try:
                radius = float(radius)
            except:
                raise ValueError('Radius must be a float or astropy.units')

        if wavelength is None:
            wavelength = np.arange(1000, 10000, 10, dtype=float) * u.Angstrom
        else:
            # TODO: Tests for whether wavelength is a 1d ndarray must be done
            pass

        self._temperature = tempeture
        self._radius = radius
        self._z = redshift
        self._input_wavelength = wavelength

        def compute_blackbody(self):
            # sed = np.zeros_like(self._input_wavelength)
            # sed = (2.0 * CGS_H * np.pi * CGS_C**2.0) / (x * 1e-8)**5
            # sed /= np.exp(CGS_H * CGS_C / (x * 1e-8 * CGS_K * T)) - 1.0
            # sed *= 1e-8
            # sed *= (R * 1e14)**2
            # sed /= (lumDis)**2
            # sed /= (1.0 + z)
            # sed *= np.interp(x, absWavelenght, absBandpass, left=0.4, right=1.0)

                return sed
