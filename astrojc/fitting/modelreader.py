from astropy.io import fits
from os import path
import numpy as np

from ..logging import log


def _float_equal(f1, f2, precision):
    if abs(float(f1) - float(f2)) < precision:
        return True
    else:
        return False


def read_castelli_model(folder, metallicity, teff, logg, catalog='catalog.fits'):
    """Read Castelli models from a folder, for a given metallicity, Teff, logg
    and a Castelli fits catalog.
    """
    c = fits.open(path.join(folder, catalog))
    fname = None
    for i in c[1].data:
        t, m, g = [np.float(k) for k in i[0].split(',')]
        if _float_equal(t, teff, 1.0) and \
           _float_equal(m, metallicity, 0.1) and _float_equal(g, logg, 0.1):
            fname = i[1]
            log.info('Loading Castelli model %s' % fname)
            pass
    if fname is None:
        log.debug('teff = %f logg=%f z=%f' % (t, g, m))
        raise ValueError("File with these parameters doesn't exists")

    f = fits.open(path.join(folder, fname.split('[')[0]))[1].data
    col = fname.split('[')[1].replace(']', '')
    return f['WAVELENGTH'], f[col]
