import numpy as np
import astropy.cosmology


class Spectrum():
    """
    Generic spectrum class

    Spectrum class for storying, loading and manipulating spectral data from
    any astronomical source. Build in operations include redshifting,
    synthesising photometry and resampling of spectra.

    Either a file name or arrays containing wavelength and flux values must be
    provided. If both get provided the file will still be read but the inputs
    will be overwritten by the values supplied.

    This class is currently expecting to receive wavelengths in Angstroms and
    fluxes in erg/s/cm^2/A. This may be generalised in the future. While the
    class will continue operating all magnitues and other statistics provided
    will be completely non sensical without the correct units.

    Parameters
    ----------
    file_name : str, optional
        Path to an ASCII file containing a spectrum. The input is expected to
        have a minimum of two columns where the first column contains
        wavelengths in angstroms and the second column fluxes in units of
        erg/s/cm^2/A.

    wavelength : array-like, optional
        Wavelength array for the specta, must be a 1d array of floats with the
        same dimention as the corresponding array of fluxes. Wavelengths must
        be provided in Angstroms. Used if an input ASCII file is not provided.

    flux : array-like, optional
        Flux array for the specta, must be a 1d array of floats with the
        same dimention as the corresponding array of wavelengths.
        Flux must be provided in the units of erg/s/cm^2/A.
        Used if an input ASCII file is not provided.

    redshift : float, optional
        Redshift of the supernova. If no value is provided z=0 will be assumed
        and luminosity_distance will be set to a default value of 10pc.

    cosmology : `astropy.cosmology.FlatLambdaCDM`, optional
        Cosmology object used to compute luminosity distances. If no value is
        provided Planck15 cosmology will be assumed by default.
    """

    def __init__(self, file_name=None, wavelength=None, flux=None,
                 redshift=0, cosmology=None):
        if file_name is not None:
            self.load_from_file(file_name)
        elif (wavelength is None) and (flux is None):
            raise IOError("""Either file_name or wavelength and
                          flux must be provided""")

        if wavelength is not None:
            try:
                iter(wavelength)
                wavelength = np.array(wavelength)
            except TypeError:
                print('wavelength must be 1d array-like')

            try:
                wavelength = wavelength.astype(float)
            except ValueError:
                print('wavelength array must numeric')

            self._wavelength = wavelength

        if flux is not None:
            try:
                iter(flux)
                flux = np.array(flux)
            except TypeError:
                print('flux must be 1d array-like')

            try:
                flux = flux.astype(float)
            except ValueError:
                print('flux array must numeric')

            self._flux = flux

        if not self._wavelength.shape == self._flux.shape:
            raise ValueError('wavelength and flux must be have same shape')

        if not self._wavelength.shape[0] == self._wavelength.size:
            raise ValueError('wavelength and flux must be 1d arrays')

        if cosmology is None:
            self.__cosmology = astropy.cosmology.Planck15
        elif type(cosmology) == astropy.cosmology.core.FlatLambdaCDM:
            self.__cosmology = cosmology
        else:
            raise TypeError("""cosmology must be of type
                            astropy.cosmology.core.FlatLambdaCDM""")

        self._z = redshift
        self._lum_distance = self.__cosmology.luminosity_distance(self._z)

    def load_from_file(self, file_name):
        """
        Load spectrum from file

        Paramaters
        ----------
        file_name : str
            Path of the ASCII spectrum input file. It must contain at least
            two columns; wavelength and flux.

        Returns
        -------
        None
        """
        try:
            arr = np.loadtxt(file_name, unpack=True)
        except IOError:
            print('Could not read {}'.format(file_name))

        if arr.shape[0] > 1:
            arr = arr[0:2]
            try:
                arr = arr.astype(float)
            except:
                raise TypeError('ASCII input file must contain floats')

        self._wavelength = arr[0]
        self._flux = arr[1]

    def adjust_redshift(self, new_redshift):
        """
        Move the spectrum to a new input redshift

        Moves the redshift of the spectrum by adjusting the wavelength and
        flux by factors of (1 + new_redshift) / (1 + old_redshift) and its
        inverse respectively. Flux is also adjusted by a square of the ratios
        of the luminosity distances at both redshifts.

        Parameters
        ----------
        new_redshift : float
            New redshift of the spectrum

        Returns
        -------
        None
        """
        redshift_shift = (1 + new_redshift) / (1 + self._z)
        self._wavelength *= redshift_shift
        self._flux /= redshift_shift
        self._z = new_redshift

        new_lum_distance = self.__cosmology.luminosity_distance(new_redshift)
        distance_factor = (self._lum_distance / new_lum_distance)**2
        self._flux /= distance_factor.value
        self._lum_distance = new_lum_distance

    def synthesis_photometry(self, filter_names, filters):
        """
        Make synthetic photometry from a spectum

        Parameters
        ----------
        filter_name : str or array-like
            Value or an array of names of filters at which the synthetic
            photometry is to be calculated.

        filters : `spek.Filters`
            Filters object that must be preloaded with the filters passed as
            filter_name

        Returns
        -------
        synthetic_flux : ndarray
            Array of synthetic fluxes matching the input filter_name
        """
        filter_names = np.array(filter_names)
        synthetic_flux = np.zeros(filter_names.size)

        for i, flt in enumerate(filter_names):
            bandpass = np.interp(self._wavelength,
                                 filters[flt].wavelength,
                                 filters[flt].bandpass, left=0, right=0)
            flux = self._flux * bandpass
            synthetic_flux[i] = (np.trapz(flux, x=self._wavelength) /
                                 filters[flt].area)

        return synthetic_flux

    def synthesis_magnitudes(self, filter_name, filters):
        """
        Calculate synthetic magnitudes for the spectrum using provided filters

        Parameters
        ----------
        filter_name : str or array-like
            Value or an array of names of filters at which the synthetic
            photometry is to be calculated.

        filters : `spek.Filters`
            Filters object that must be preloaded with the filters passed as
            filter_name

        Returns
        -------
        synthetic_mag: ndarray
            Array of synthetic magnitudes matching the input filter_name
        """
        synthetic_flux = self.synthesis_photometry(filter_name, filters)
        synthetic_mag = np.zeros_like(synthetic_flux)

        for i, flt in enumerate(filter_name):
            synthetic_mag[i] = (-2.5*np.log10(synthetic_flux[i]) -
                                filters[flt].ab_zero_point)

        return synthetic_mag
