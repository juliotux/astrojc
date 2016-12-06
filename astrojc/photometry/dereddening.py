'''Dereddening of a SED or spectrum, using literature reddeing curves.'''

import numpy as np
from numpy.polynomial.polynomial import polyval
from scipy.interpolate import splrep, splev

__all__ = ['unred_ccm89', 'unred_gcmlw03']

try:
    from numba import vectorize
    numbajit = True
    target = 'cpu'
except:
    numbajit = False

# Gordon et al. 2003, Apj, 594, 279
_gordon_ref = np.array(list(zip([0.455, 0.606, 0.8, 1.235, 1.538, 1.818, 2.273, 2.703,
                          3.375, 3.625, 3.875, 4.125, 4.375, 4.625, 4.875,
                          5.125, 5.375, 5.625, 5.875, 6.125, 6.375, 6.625,
                          6.875, 7.125, 7.375, 7.625, 7.875, 8.125, 8.375, 8.625],
                         [0.016, 0.169, 0.131, 0.567, 0.801, 1.0, 1.374, 1.672,
                          2.0, 2.22, 2.428, 2.661, 2.947, 3.161, 3.293, 3.489,
                          3.637, 3.866, 4.013, 4.243, 4.472, 4.776, 5.0, 5.272,
                          5.575, 5.795, 6.074, 6.297, 6.436, 6.992],
                         [0.101, 0.097, 0.299, np.nan, np.nan, 1.0, 1.349, 1.665,
                          1.899, 2.067, 2.249, 2.447, 2.777, 2.922, 2.921, 2.812,
                          2.805, 2.863, 2.932, 3.06, 3.11, 3.299, 3.408, 3.515,
                          3.67, 3.862, 3.937, 4.055, 3.969, np.nan],
                         [0.030, 0.186, 0.257, np.nan, np.nan, 1.0, 1.293, 1.518,
                          1.786, 1.969, 2.149, 2.391, 2.771, 2.967, 2.846, 2.646,
                          2.565, 2.566, 2.598, 2.607, 2.668, 2.787, 2.874, 2.983,
                          3.118, 3.231, 3.374, 3.366, 3.467, np.nan])),
                         dtype=np.dtype([('x', 'f8'), ('SMCBAR', 'f8'),
                                         ('LMC2', 'f8'), ('AVGLMC', 'f8')]))

def unred_gcmlw03(wavelength, flux, loc, Av, Rv=3.1, mag=False):
    '''
    Deredden a flux vector using the reddening curves from Gordon et al. 2003,
    Apj, 594, 279.

    Parameters:
    -----------
        wavelength : array_like of float
            The wavelength to calculate, in microns.
        flux : array_like or float
            The flux data to be derreded.
        loc : (`SMCBAR`, `LMC2`, `AVGLMC`)
            The extiction law (location) of the curve to be used to unredden the
            data.
        Av : float
            The reference Av.
        Rv : float, optional
            The ratio Av/E(B-V). For non anomalous extinction law, it is 3.1.
            But the following values are commonly assumed:
                - MW Difuse          rv = 3.1
                - MW Dense           rv = 5.0
                - LMC Average        rv = 3.41
                - LMC2Supershell     rv = 2.76
                - SMC Bar            rv = 2.74
        mag : bool, optional
            If true, the flux data will be interpreted in units of magnitude. If
            false, will be interpreted as a linear scale.
    Retruns:
    --------
        dered_flux : array_like
            The derreded flux data.
        AlAv : array_like
            The reddening array calculated for the data, relative to Av.
    '''
    x = 1/np.array(wavelength)
    locs = ['SMCBAR', 'LMC2', 'AVGLMC']
    if loc not in locs:
        raise ValueError('Cannot understood {} extinction law.'
                         'The values are: {}'.format(loc, locs))
    x = 1/np.array(wavelength)
    wr, ar = _gordon_ref['x'], _gordon_ref[loc]
    filt = np.isnan(ar)
    wr, ar = wr[np.where(~filt)], ar[np.where(~filt)]
    spline = splrep(wr, ar, k=3)
    Al = Av*splev(x, spline)
    if mag:
        dered_flux = np.array(flux) - Al
    else:
        dered_flux = np.array(flux)*10.0**(0.4*Al)
    return dered_flux, Al

# Cardelli et al. 1989, Apj, 345, 245
def unred_ccm89(wavelength, flux, Av, Rv=3.1, mag=False):
    '''
    Deredden a flux vector using the parameterization from Cardelli, Clayton,
    and Mathis (1989 ApJ. 345, 245). Coefficients from from O'Donnell (1994).

    Adapted from the IDL routine ccm_unred.pro.

    Parameters:
    -----------
        wavelength : array_like of float
            The wavelength to calculate, in microns.
        flux : array_like or float
            The flux data to be derreded.
        Av : float
            The reference Av.
        Rv : float
            The ratio Av/E(B-V). For non anomalous extinction law, it is 3.1.
            But the following values are commonly assumed:
                - MW Difuse          rv = 3.1
                - MW Dense           rv = 5.0
                - LMC Average        rv = 3.41
                - LMC2Supershell     rv = 2.76
                - SMC Bar            rv = 2.74
        mag : bool, optional
            If true, the flux data will be interpreted in units of magnitude. If
            false, will be interpreted as a linear scale.
    Retruns:
    --------
        dered_flux : array_like
            The derreded flux data.
        Al : array_like
            The reddening array calculated for the data
    '''
    x = 1/np.array(wavelength)
    ax = np.zeros(len(x))
    bx = np.zeros(len(x))
    ax.fill(np.nan)
    bx.fill(np.nan)

    #IR
    filt = np.where((x >= 0.3) & (x < 1.1))
    if len(filt) > 0:
        ax[filt] = 0.574 * x[filt]**(1.61)
        bx[filt] = -0.527 * x[filt]**(1.61)
    #NIR/Optical
    filt = np.where((x >= 1.1) & (x < 3.3))
    if len(filt) > 0:
        ca = [1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505]
        cb = [0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347]
        ax[filt] = polyval(x[filt], ca)
        bx[filt] = polyval(x[filt], cb)
    #UV
    filt = np.where((x >= 3.3) & (x < 8.0))
    if len(filt) > 0:
        y = x[filt]
        fa = np.zeros(len(filt))
        fb = np.zeros(len(filt))
        filt1 = np.where(y >= 5.9)
        if len(filt1) > 0:
            caf = [0, 0, -0.04473, -0.009779]
            cbf = [0, 0, 0.2130, 0.1207]
            fa[filt1] = polyval(y - 5.9, caf)
            fb[filt1] = polyval(y - 5.9, cbf)
        ax[filt] = 1.752 - 0.316*y - (0.104/((y-4.67)**2 + 0.341)) + fa
        bx[filt] = -3.090 + 1.825*y + (1.206/((y-4.62)**2 + 0.263)) + fb
    #FUV
    filt = np.where((x >= 8.0) & (x < 11.0))
    if len(filt) > 0:
        ca = [-1.073, -0.628, 0.137, -0.070]
        cb = [13.670, 4.257, -0.420, 0.374]
        ax[filt] = polyval(x[filt], ca)
        bx[filt] = polyval(x[filt], cb)

    Al = Av*(ax + bx/Rv)
    if mag:
        dered_flux = np.array(flux) - Al
    else:
        dered_flux = np.array(flux)*10.0**(0.4*Al)
    return dered_flux, Al
