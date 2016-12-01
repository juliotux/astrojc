from astropy.io import fits
import numpy as np

__all__ = ['iraf_fits_loader', 'xshooter_sflx_loader']

def xshooter_sflx_loader(fname):
    '''
    Reads an XSHOOTER sflx file.
    '''
    f = fits.open(fname)
    return f[1].data['wave'][0], f[1].data['flux'][0], f[0].header

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

    return wavelength, flux, header

