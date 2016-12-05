import numpy as np

__all__ = ['rv_correct', 'Spectrum']

def rv_correct(wave, rv):
    '''
    Correct a wavelenght array with radial velocity, in km/s

    Parameters:
    -----------
    wave : array_like
        Array of wavelength to be corrected by radial velocity.
    rv : float
        Radial velocity of the object.

    Returns:
    --------
    array_like : The radial velocity corrected wavelength array.
    '''
    return wave*(1 - rv/3e5)

class Spectrum(object):
    '''Generic class for 1D spectrum.'''
    #TODO: Implement a way to calculate wcs. Raise error if WCS problem (like ESO tfits)
    #TODO: Implement quantities
    def __init__(self, wave, flux, header, rv=0):
        self._wave = wave
        self._flux = flux
        self._header = header
        self._rv = rv

    @property
    def wave(self):
        return self.fix_rv(self._rv)

    @wave.setter
    def wave(self, wave):
        self._wave = wave

    @property
    def flux(self):
        return self._flux

    @flux.setter
    def flux(self, flux):
        self._flux = flux

    @property
    def rv(self):
        return self._rv

    @rv.setter
    def rv(self, rv):
        self._rv = rv

    @property
    def wcs(self):
        raise NotImplementedError

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, header):
        self._header = header

    def fix_rv(self, rv=None):
        if rv is None:
            rv = self._rv
        return rv_correct(self._wave, rv)

    def append(self, spectrum):
        self._wave = np.append(self._wave, spectrum.wave)
        self._flux = np.append(self._flux, spectrum.flux)
