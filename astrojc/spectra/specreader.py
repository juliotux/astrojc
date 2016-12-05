from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import numpy as np

from .spectrum import Spectrum

__all__ = ['iraf_fits_loader', 'xshooter_sflx_loader', 'xshooter_full_loader']

def xshooter_sflx_loader(fname):
    '''
    Reads an XSHOOTER sflx file.
    '''
    #TODO: when Spectrum unit implemented, implement unit here.
    f = fits.open(fname)
    return Spectrum(f[1].data['wave'][0],# * u.Unit(f[1].header['TUNIT1']),
                    f[1].data['flux'][0],# * u.Unit(f[1].header['TUNIT2']),
                    f[0].header)

def xshooter_full_loader(blue_fname = None, vis_fname = None, ir_fname = None,
                         blue_correct_factor = 1.0,
                         vis_correct_factor = 1.0,
                         ir_correct_factor = 1.0):
    '''
    Reads spec from all the 3 arms of XSHOOTER and return one single spectra.
    The flux difference of arms can be corrected using the correct_factor for
    each one, that is a number that will multiply the flux of the arm.
    '''
    #TODO: implement a header joining
    spec = Spectrum(np.zeros(0), np.zeros(0), header=None)

    for s,cf,ll,lh in zip([blue_fname, vis_fname, ir_fname],
                          [blue_correct_factor, vis_correct_factor, ir_correct_factor],
                          [ 310.0,  555.6,  1020.0],
                          [ 555.6, 1020.0,  2480.0]):
        if s is not None:
            sp = xshooter_sflx_loader(s)
            #TODO: when Spectrum unit implemented, implement unit here.
            filt = (sp.wave >= ll) & (sp.wave < lh)
            sp._wave = sp._wave[np.where(filt)]
            sp._flux = cf*sp._flux[np.where(filt)]
            spec.append(sp)
    return spec

def iraf_fits_loader(fname):
    '''
    Reads standart IRAF spectrum.
    '''
    sp = fits.open(fname)
    header = sp[0].header

    wcs = WCS(header)
    #make index array
    index = np.arange(header['NAXIS1'])

    wavelength = wcs.wcs_pix2world(index[:,np.newaxis], 0)
    wavelength = wavelength.flatten()
    flux = sp[0].data

    return Spectrum(wavelength, flux, header)

